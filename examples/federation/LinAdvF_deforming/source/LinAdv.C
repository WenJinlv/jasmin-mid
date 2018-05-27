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
#include "NodeData.h"
#include "SideData.h"
#include "FederalDeformingPatchGeometry.h"
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
               tbox::Pointer<tbox::Database> input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
#endif

   d_object_name = object_name;

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
  
   // 当前值上下文, 新值上下文, 演算上下文.
   d_current = variable_db->getContext("CURRENT");
   d_new     = variable_db->getContext("NEW");
   d_scratch = variable_db->getContext("SCRATCH");

   // 定义变量.
   d_uval = variable_db->getVariable("uval");
   if(d_uval.isNull()) d_uval = new pdat::CellVariable<NDIM,double>("uval",1);
   d_uval_current_id = 
        variable_db->registerVariableAndContext(d_uval, d_current);
   d_uval_new_id = 
        variable_db->registerVariableAndContext(d_uval, d_new );
   d_uval_scratch_id = 
        variable_db->registerVariableAndContext(d_uval, d_scratch, 1);
    
   d_flux = variable_db->getVariable("flux");
   if(d_flux.isNull()) d_flux = new pdat::SideVariable<NDIM,double>("flux",1);
   d_flux_new_id = variable_db->registerVariableAndContext(d_flux, d_new );

   d_coords = variable_db->getVariable("coords");
   if(d_coords.isNull()) d_coords = new pdat::NodeVariable<NDIM,double>("coords",NDIM);
   d_coords_id = 
        variable_db->registerVariableAndContext(d_coords, d_current);
}

/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void LinAdv::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string &intc_name = intc->getName();

   if ( intc_name=="ALLOC_NEW_DATA") {  // 内存构件.
         intc->registerPatchData(d_uval_new_id);   
         intc->registerPatchData(d_flux_new_id);  
    }else if  ( intc_name=="ALLOC_SCRATCH_DATA") {  // 内存构件.
         intc->registerPatchData(d_uval_scratch_id);
    } else if ( intc_name=="INIT_SET_VALUE") {  // 初值构件. 
         intc->registerInitPatchData(d_uval_current_id);
         intc->registerInitPatchData(d_coords_id);

    } else if(  intc_name=="STEP_SIZE" ) {  // 步长构件.

    } else if ( intc_name=="COMPUTE_FLUX") {  // 数值构件: 计算通量，需要u的影像区数据.
          intc->registerRefinePatchData(d_uval_scratch_id,
                                        d_uval_current_id);

    } else if ( intc_name=="CONSER_DIFF") {  // 数值构件: 计算守恒量, 不需要影像区数据.

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
      tbox::Pointer< pdat::NodeData<NDIM,double> > coords =
                                    patch.getPatchData(d_coords_id);

      const hier::IntVector<NDIM> &ghost_cells = uval_current->getGhostCellWidth();

      const tbox::Pointer<geom::FederalDeformingPatchGeometry<NDIM> > 
                  pgeom = patch.getPatchGeometry();

      const int fn = pgeom->getFederalNumber();

      const hier::Index<NDIM> &ifirst=patch.getBox().lower();
      const hier::Index<NDIM> &ilast =patch.getBox().upper();

      if (fn==1) {
/*
      initsetvalue_(
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
*/
      initcoords1_(
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
                    coords->getPointer());

         uval_current->fillAll(1.0);
      }else if (fn==0) {
         uval_current->fillAll(0.0);

      initcoords0_(
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
                    coords->getPointer());

      }
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

   return(1.0);

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
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="COMPUTE_FLUX" || intc_name=="CONSER_DIFF" );
#endif

   const tbox::Pointer<geom::FederalDeformingPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();

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
         advancepatch_(d_x_velocity, dt,
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

         const int fn = pgeom->getFederalNumber();

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
	 conspatch_(
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
      if (fn==0 && !initial_time) {
         uval_current->fillAll(2.0);
      }

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
   assert(intc_name=="COMPUTE_FLUX" || intc_name == "INIT_SET_VALUE"); 
#endif

   NULL_USE(fill_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > uval 
                          = patch.getPatchData(d_uval_scratch_id);

   const hier::IntVector<NDIM> &ghost_cells = uval->getGhostCellWidth();

   const tbox::Pointer< hier::FederalPatchGeometry<NDIM> > pgeom =
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
 *  从输入数据库读入数据. 
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

   javis_writer->registerNodeCoordinates(d_coords_id);
   javis_writer-> registerPlotQuantity("UUU", "SCALAR", d_uval_current_id);

}


