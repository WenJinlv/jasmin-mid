//
// 文件名: LinAdv.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 $
// 描述  : 线性对流问题的网格片时间积分算法类的实现.
//

#include "LinAdv.h"
#include "LinAdvFort.h"

#include "CellData.h"
#include "SideData.h"
#include "UniRectangularPatchGeometry.h"
#include "tbox/RestartManager.h"
using namespace JASMIN;

#include <assert.h>
using namespace std;


/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
LinAdv::LinAdv(const string& object_name,
               tbox::Pointer<tbox::Database> input_db,
               tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
#endif

   d_object_name = object_name;
   d_grid_geometry = grid_geom;

   // 读取从输入文件或重启动文件读入数据.
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   } else {
      getFromInput(input_db);
   }

   // 注册变量和数据片.
   registerModelVariables();

   // 注册为重启动对象.
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
LinAdv::~LinAdv()
{
tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};


/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void LinAdv::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   d_uval    = new pdat::CellVariable<NDIM,double>("uval",1);  // 中心量，深度为1.
   d_flux    = new pdat::SideVariable<NDIM,double>("flux",1);  // 边心量，深度为1.
  
   // 当前值上下文, 新值上下文, 演算上下文.
   d_current = variable_db->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = variable_db->getContext("SCRATCH");

   // 数据片<uval,current>: 存储当前时刻的 u 值, 影像区宽度为0.
   d_uval_current_id = variable_db->registerVariableAndContext(d_uval, d_current);

   // 数据片<uval,new>: 存储新时刻的 u 值, 影像区宽度为0.
   d_uval_new_id = variable_db->registerVariableAndContext(d_uval, d_new );

   // 数据片<flux,new>: 存储新时刻的 f 值, 影像区宽度为0.
   d_flux_new_id = variable_db->registerVariableAndContext(d_flux, d_new );

   // 数据片<uval,scratch>: 存储变量 u 的演算值, 影像区宽度为1.
   d_uval_scratch_id = variable_db->registerVariableAndContext(d_uval, d_scratch, 1);

}


/*************************************************************************
 *                                                                        
 * 注册绘图量.
 *                                                                  
 *************************************************************************/
void LinAdv::registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(javis_writer.isNull()));
#endif

   javis_writer-> registerPlotQuantity("UUU", "SCALAR", d_uval_current_id);

}


/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void LinAdv::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{

   const string &intc_name = intc->getName();

   if ( intc_name=="ALLOC_NEW_DATA") {  // 内存构件.
         intc->registerPatchData(d_uval_new_id);   
         intc->registerPatchData(d_flux_new_id);  

    } else if ( intc_name=="ALLOC_SCRATCH_DATA") {  // 内存构件.
         intc->registerPatchData(d_uval_scratch_id);  

    } else if ( intc_name=="INIT_SET_VALUE") {  // 初值构件. 
         intc->registerInitPatchData(d_uval_current_id,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                     d_uval_scratch_id);

    } else if(  intc_name=="STEP_SIZE" ) {  // 步长构件.

    } else if ( intc_name=="COMPUTE_FLUX") {  // 数值构件: 计算通量，需要u的影像区数据.
         intc->registerRefinePatchData(d_uval_scratch_id, // 目的数据片
                                       d_uval_current_id, // 同层源数据片
                                      "CONSERVATIVE_LINEAR_REFINE", // 细化插值算子
                                       d_uval_scratch_id, // 演算数据片
                                       d_uval_current_id, // 粗层源数据片1（时间插值时需要）
                                       d_uval_new_id);    // 粗层源数据片2（时间插值时需要）

    } else if ( intc_name=="CONSER_DIFF") {  // 数值构件: 计算守恒量, 不需要影像区数据.

    } else if ( intc_name=="TAG_CELLS_GEOM") { // 数值构件: 标记待重建的计算区域.

    } else if ( intc_name=="TAG_CELLS") {    // 数值构件: 标记细化单元, 需要u的影像区数据.
          intc->registerRefinePatchData(d_uval_scratch_id,
                                        d_uval_current_id,
                                       "CONSERVATIVE_LINEAR_REFINE");  

    } else if ( intc_name=="SYNC_INTERFACE") {  // 同步构件: 粗化通量.
           intc->registerCoarsenPatchData(d_flux_new_id,
                                         "CONSERVATIVE_COARSEN");

    } else if ( intc_name=="SYNC_OVERLAY") {    // 同步构件: 粗化守恒量.
           intc->registerCoarsenPatchData(d_uval_new_id,
                                         "CONSERVATIVE_COARSEN",
                                          d_uval_new_id);

    } else if ( intc_name=="COPY_SOLUTION") {  // 复制构件.
          intc->registerCopyPatchData(d_uval_current_id,
                                      d_uval_new_id);   

    } else {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
    }

}

/*************************************************************************
 *
 *  初始化数据片（支持初值构件）.
 *
 ************************************************************************/
void LinAdv::initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT_SET_VALUE");
#endif

   (void) time;

   if (initial_time) {
      tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
                                    patch.getPatchData(d_uval_current_id);

      const hier::IntVector<NDIM> &ghost_cells = uval_current->getGhostCellWidth();

      const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const hier::Index<NDIM> &ifirst=patch.getBox().lower();
      const hier::Index<NDIM> &ilast =patch.getBox().upper();

      initsetvalue_(dx,xlo,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    uval_current->getPointer());
   }

}

/*************************************************************************
 *
 *  计算稳定性时间步长（支持步长构件）.
 *
 ************************************************************************/
double LinAdv::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="STEP_SIZE");
#endif

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom
                       = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();

   return 0.8*dx[0]/d_x_velocity;

}

/*************************************************************************
 *
 * 完成单个网格片上的数值计算（支持数值构件）.
 *
 * 该函数基于显式迎风格式，实现通量和守恒量的计算.
 ************************************************************************/
void LinAdv::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* xlo = pgeom->getXLower();

   const hier::Index<NDIM> &ifirst=patch.getBox().lower();
   const hier::Index<NDIM> &ilast =patch.getBox().upper();

   if(intc_name=="COMPUTE_FLUX" ) {  

         // 计算通量

         tbox::Pointer< pdat::SideData<NDIM,double> > flux_new =
                                    patch.getPatchData(d_flux_new_id);
         tbox::Pointer< pdat::CellData<NDIM,double> > uval_scratch =
                                    patch.getPatchData(d_uval_scratch_id);

         hier::IntVector<NDIM> ghost_cells = uval_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
         // 检查通量的影像区宽度是否与FORTRAN函数(advancepatch_)中的声明一致.
	 assert(flux_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
         advancepatch_(d_x_velocity, dt, dx,xlo,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    uval_scratch->getPointer(),
                    flux_new->getPointer(0),
                    flux_new->getPointer(1)
#if (NDIM==3)
                    ,flux_new->getPointer(2)
#endif
                    );

   } else if(intc_name=="CONSER_DIFF" ) {

         // 计算守恒量
         
         tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
                                    patch.getPatchData(d_uval_current_id);
	 tbox::Pointer< pdat::CellData<NDIM,double> > uval_new =
                                    patch.getPatchData(d_uval_new_id);          
	 tbox::Pointer< pdat::SideData<NDIM,double> > flux_new =
                                    patch.getPatchData(d_flux_new_id);

	 const hier::IntVector<NDIM> &ghost_cells =
                                    uval_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
         // 检查通量的影像区宽度是否与FORTRAN函数(conspatch_)中的声明一致.
	 assert(flux_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
	 conspatch_(dx, xlo,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    uval_current->getPointer(),
                    flux_new->getPointer(0),
                    flux_new->getPointer(1),
#if (NDIM==3)
                    flux_new->getPointer(2),
#endif
                    uval_new->getPointer() );

   } else if(intc_name=="TAG_CELLS" ) {

         tbox::Pointer< pdat::CellData<NDIM,double> > uval_scratch =
                                    patch.getPatchData(d_uval_scratch_id);
	 tbox::Pointer< pdat::CellData<NDIM,int> > tag =
                                    patch.getPatchData(this->getTagIndex());

	 const hier::IntVector<NDIM> &ghost_cells = uval_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
	 // 检查通量的影像区宽度是否与FORTRAN函数(conspatch_)中的声明一致.
	 assert(tag->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif

	 tagpatch_(dx,xlo,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    uval_scratch->getPointer(),
                    tag->getPointer()
                       );

   } else if(intc_name=="TAG_CELLS_GEOM" ) {
	 tbox::Pointer< pdat::CellData<NDIM,int> > tag =
                                    patch.getPatchData(this->getTagIndex());

         const double* xlo_global = d_grid_geometry->getXLower();
         const double* xup_global = d_grid_geometry->getXUpper();
         
	 const hier::IntVector<NDIM> &ghost_cells = tag->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
	 // 检查通量的影像区宽度是否与FORTRAN函数(conspatch_)中的声明一致.
	 assert(tag->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif

	 tagpatchforgeom_(xlo_global, xup_global,
                    dx,xlo,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    tag->getPointer()
                       );

   } else { 
           TBOX_ERROR("\n::computeOnPatch() : component "     
                           << intc_name <<" is not matched. "<<endl);
   }
}


/*************************************************************************
 *
 *  填充物理边界条件.
 *
 ************************************************************************/
void LinAdv::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="COMPUTE_FLUX" || intc_name == "INIT_SET_VALUE" ||
          intc_name=="TAG_CELLS" );
#endif

   NULL_USE(fill_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > uval 
                          = patch.getPatchData(d_uval_scratch_id);

   const hier::IntVector<NDIM> &ghost_cells = uval->getGhostCellWidth();

   const tbox::Pointer< hier::PatchGeometry<NDIM> > pgeom =
                                       patch.getPatchGeometry();

   const hier::Box<NDIM> &patch_box = patch.getBox();

   for(int btype=1; btype<=NDIM; btype++) {

      // 获取物理边界影像区
      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,
                                  loc_indexes,
                                  btype, 
                                  patch_box,
                                  ghost_cells);

      // 填充物理边界影像区.
#if ( NDIM == 2 )  // 二维情形
      if( btype == EDGE2D_BDRY_TYPE ) {
          for (int num=0; num < fill_boxes.size(); num++) {
                if ( loc_indexes[num] == XLO ) {
                     uval->fillAll(2.0,fill_boxes(num));
                } else {
                     uval->fillAll(0.0,fill_boxes(num));
                }
          }
      } else {
          for (int num=0; num < fill_boxes.size(); num++) {
                if ( loc_indexes[num] == XLO_YLO 
                     || loc_indexes[num] == XLO_YHI ) {
                     uval->fillAll(2.0,fill_boxes(num));
                } else {
                     uval->fillAll(0.0,fill_boxes(num));
                }
	  }
      }
#elif ( NDIM == 3 )  // 三维情形
      if( btype == FACE3D_BDRY_TYPE ) {
          for (int num=0; num < fill_boxes.size(); num++) {
                if ( loc_indexes[num] == XLO ) {
                     uval->fillAll(2.0,fill_boxes(num));
                } else {
                     uval->fillAll(0.0,fill_boxes(num));
                }
          }
      } else if( btype == EDGE3D_BDRY_TYPE ) {
          for (int num=0; num < fill_boxes.size(); num++) {
                if ( loc_indexes[num] == XLO_YLO 
                     || loc_indexes[num] == XLO_YHI ) {
                     uval->fillAll(2.0,fill_boxes(num));
                } else {
                     uval->fillAll(0.0,fill_boxes(num));
                }
	  }
      } else {
          for (int num=0; num < fill_boxes.size(); num++) {
                if ( loc_indexes[num] == XLO_YLO_ZLO
                     || loc_indexes[num] == XLO_YLO_ZHI
                     || loc_indexes[num] == XLO_YHI_ZLO
                     || loc_indexes[num] == XLO_YHI_ZHI ) {
                     uval->fillAll(2.0,fill_boxes(num));
                } else {
                     uval->fillAll(0.0,fill_boxes(num));
                }
	  }
      }
#endif
   }

}

/*************************************************************************
 *
 *  从输入数据库读取数据. 
 *
 ************************************************************************/
void LinAdv::getFromInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("constant_x_velocity")) {
       d_x_velocity = db->getDouble("constant_x_velocity");
       TBOX_ASSERT( d_x_velocity > 0.0 );
   } else {
       TBOX_ERROR(d_object_name << ": "
              << " No key `constant_x_velocity' found in data." << endl);
   }
}

/*************************************************************************
 *
 *  输出数据成员到重启动数据库.
 *
 ************************************************************************/
void LinAdv::putToDatabase(tbox::Pointer<tbox::Database> db)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putDouble("d_x_velocity", d_x_velocity);

}


/*************************************************************************
 *
 *  从重启动数据库读取数据. 
 *
 ************************************************************************/
void LinAdv::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   d_x_velocity = db->getDouble("d_x_velocity");

}

