//
// 文件名:  FDMBSingDemo.C
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 170 $
// 修改  :  $Date: 2007-06-27 08:44:22  $
// 描述  :  单个网格片上的数值计算子程序.
//

#ifndef included_FDMBSingDemo_C
#define included_FDMBSingDemo_C

#include "FDMBSingDemo.h"

using namespace std;

#include "FederalBlockPatchGeometry.h"

#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/IEEE.h"
#include "tbox/Array.h"

//声明Fortran数值子程序的头文件.
#include "FDMBSingDemoFort.h"

//解向量的分量个数
#define NEQU  (NDIM + 2) //NDIM个速度分量+ 压力 + 密度

/*
*************************************************************************
*                                                                       *
*  FDMBSingDemo类的构造函数, 执行以下操作:                                     *
*  (1) 初始化数据成员                                                   *
*  (2) 创建定义物理模型控制方程中的数值解的变量.                        *
*  (3) 如果进行断点续算, 那么调用getFromRestart().                      *
*  (4) 调用getFromInput(), 从指定数据库中读取数据.                      *
*      该操作可能覆盖从重启动文件中读入的数据.                          *
*                                                                       *
*************************************************************************
*/
FDMBSingDemo::FDMBSingDemo(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridConstructionUtilities<NDIM,double> > grid_tool, 
         tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!grid_geom.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometry = grid_geom;
   d_grid_tool1     = grid_tool;

   /*
    * 从指定的输入数据库/重启动数据库初始化数据成员.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   // 为影像区的宽度赋值.
   d_zeroghosts = hier::IntVector<NDIM>(0);
   d_oneghosts  = hier::IntVector<NDIM>(1);
   d_twoghosts  = hier::IntVector<NDIM>(2);

   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();

   // 创建数学运算.
   d_math_ops = new math::PatchCellDataBasicOps<NDIM,double>();

}

FDMBSingDemo::FDMBSingDemo(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
         tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!grid_geom.isNull());
   TBOX_ASSERT(!grid_tool.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometry = grid_geom;
   d_grid_tool2     = grid_tool;

   /*
 *     * 从指定的输入数据库/重启动数据库初始化数据成员.
 *         */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   // 为影像区的宽度赋值.
   d_zeroghosts = hier::IntVector<NDIM>(0);
   d_oneghosts  = hier::IntVector<NDIM>(1);
   d_twoghosts  = hier::IntVector<NDIM>(2);
   
   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();
   
   // 创建数学运算.
   d_math_ops = new math::PatchCellDataBasicOps<NDIM,double>();
   
}

/*
*************************************************************************
*    FDMBSingDemo 类的空析构函数.	 	                                *
*************************************************************************
*/

FDMBSingDemo::~FDMBSingDemo() 
{ 

tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);

if(d_math_ops) delete d_math_ops;

} 

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/

void FDMBSingDemo::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");

   // 存储当前时刻的变量值: (coords,current), 影像区宽度为0.
   d_coords = variable_db->getVariable("coordinates");
   if(d_coords.isNull()) 
       d_coords = new pdat::NodeVariable<NDIM,double>("coordinates",NDIM);
   d_coords_current_id =
      variable_db->registerVariableAndContext(d_coords, d_current);
   d_coords_scratch_id =
      variable_db->registerVariableAndContext(d_coords, d_scratch, d_twoghosts);

   // 存储当前时刻的变量值: (center_var,current), 影像区宽度为0.
   d_center_var = variable_db->getVariable("center_var");
   if(d_center_var.isNull()) 
      d_center_var = new pdat::CellVariable<NDIM,double>("center_var",1);
   d_center_var_current_id =
      variable_db->registerVariableAndContext(d_center_var, d_current);
   d_center_var_new_id =
      variable_db->registerVariableAndContext(d_center_var, d_new);
   d_center_var_scratch_id =
      variable_db->registerVariableAndContext(d_center_var,
                                              d_scratch,
                                              d_oneghosts);

}

/*
*************************************************************************
*                                                                       *
* 注册JaVis数据输出器, 以采用JaVis工具对绘图文件进行后处理.           *
*                                                                       *
*************************************************************************
*/
void FDMBSingDemo::registerJaVisDataWriter(
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(javis_writer.isNull()));
#endif
   d_javis_writer = javis_writer;

   // 注册可视化量.
   #ifdef HAVE_HDF5

   if (!(d_javis_writer.isNull())) {

      d_javis_writer->registerNodeCoordinates(d_coords_current_id);

      d_javis_writer->registerPlotQuantity(
                              "center_var",
                              "SCALAR",
                              d_center_var_current_id,
                              0,1);

   }
   #endif

}
	 
/********************************************************************
*  初始化积分构件.
*********************************************************************/
void FDMBSingDemo::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件 : 为网格坐标和数值解赋初值.
        intc->registerInitPatchData(d_coords_current_id);
        intc->registerInitPatchData(d_center_var_current_id);
    } else if  ( intc_name=="ALLOC_SCRATCH_DATA") {  // 内存构件.
        intc->registerPatchData(d_coords_scratch_id);
        intc->registerPatchData(d_center_var_scratch_id);
   }else if(intc_name=="TIME_STEP_SIZE") { // 步长构件：求步长.
   }else if(intc_name=="SUMMING") { //数值构件: 求中心量的环绕累加和.
        intc->registerRefinePatchData(d_center_var_scratch_id,
                                      d_center_var_current_id);
        intc->registerRefinePatchData(d_coords_scratch_id,
                                      d_coords_current_id);
   }else if(intc_name=="ACCEPT_SOLUTION") { //复制构件: 接收数值解.
        intc->registerCopyPatchData(d_center_var_current_id,
                                    d_center_var_new_id);
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 内存构件: 为新值数据片调度内存空间.
        intc->registerPatchData(d_center_var_new_id);
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 内存构件: 为演算数据片调度内存空间.
        intc->registerPatchData(d_center_var_scratch_id);
   }else{
        TBOX_ERROR("\n::initializeComponent() : component "
                   << intc_name <<" is not matched. "<<endl);
   }    

}

/*
*************************************************************************
* 初值构件: 网格结点坐标和数值解.
*************************************************************************
*/
void FDMBSingDemo::initializePatchData( hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const bool    initial_time,
                                 const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="INIT");
   #endif

   (void) time;

   if (initial_time) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > center_var_current = 
                               patch.getPatchData(d_center_var_current_id);


#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!coords_current.isNull());
      TBOX_ASSERT(!center_var_current.isNull());
#endif

      hier::IntVector<NDIM> ghost_cells = center_var_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ghost_cells == d_zeroghosts);
#endif

      // 生成初始网格.
      if (!d_grid_tool1.isNull()) {
         d_grid_tool1->generateDeformingMeshForDomain(patch,
                                                     d_coords_current_id);
      }else if (!d_grid_tool2.isNull()) {
         d_grid_tool2->generateDeformingMeshForDomain(patch,
                                                     d_coords_current_id);
      }

      // 赋初值.
      center_var_current->fillAll(1.0);

   }

}

/*
*************************************************************************
* 步长构件: 计算并返回网格片上的稳定时间步长.  
*************************************************************************
*/
double FDMBSingDemo::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="TIME_STEP_SIZE");
#endif

   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(initial_time);
   NULL_USE(flag_last_dt);
   NULL_USE(last_dt);
   NULL_USE(intc_name);

   return(1.0);

}

/********************************************************************
*  实现数值构件.
*********************************************************************/
void FDMBSingDemo::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="SUMMING");
   #endif

   if(intc_name=="SUMMING") {
      tbox::pout << "\nintc=" << intc_name <<",computeOnPatch:" << endl;
      summing_density(patch,time,dt,initial_time);
   }

}

// 求密度的和.
void FDMBSingDemo::summing_density(hier::Patch<NDIM>& patch,
                                     const double  time,
                                     const double  dt,
                                     const bool    initial_time)
{
   NULL_USE(time);
   NULL_USE(dt);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch =
                               patch.getPatchData(d_center_var_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_new =
                               patch.getPatchData(d_center_var_new_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!center_var_scratch.isNull());
   TBOX_ASSERT(!center_var_new.isNull());
#endif

   hier::IntVector<NDIM> cghost = center_var_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(cghost == d_oneghosts);
#endif

   const tbox::Pointer< hier::FederalBlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();

   hier::Box<NDIM> patch_box = patch.getBox();

   tbox::pout<< "federal_" << pgeom->getFederalNumber()
             << ", block_"   << pgeom->getBlockNumber()
             << ", patch_box = " << patch_box << endl;

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getFederalRegularBoundaryFillBoxes(fill_boxes,loc_indexes,
                                       btype, patch_box, cghost);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);
 
         tbox::pout << "\tfederal_regular_fill_box="<< fill_box 
                    << ", type=" << btype
                    << ", loc=" << loc_indexes[num]
                    << endl;

      }

   }

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getBlockRegularBoundaryFillBoxes(fill_boxes,loc_indexes,
                                       btype, patch_box, cghost);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);
 
         tbox::pout << "\tblock_regular_fill_box="<< fill_box << ", type=" << btype
                    << ", loc=" << loc_indexes[num]
                    << endl;

      }

   }
   
   for(int btype=2; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      tbox::Array< tbox::Array< tbox::Pointer <hier::Patch<NDIM> > > > patches;
      tbox::Array< tbox::Array< int > > blocks;
      pgeom->getBlockInternalSingularBoundaryFillBoxes(fill_boxes,
                                               loc_indexes,
                                               patches,
                                               blocks,
                                               btype,
                                               patch_box, 
                                               cghost);

      const int myid = pgeom->getBlockNumber();

      for (int num=0; num < loc_indexes.size(); num++) {

         center_var_scratch->fillAll(0.0,fill_boxes(num));

         tbox::pout << "\tblock_internal_fill_box="<< fill_boxes(num)
                    << ", type=" << btype
                    << ", loc=" << loc_indexes[num]
                    << endl;
        
         const int nblocks= blocks[num].size();

         int myloc = -1;
         for(int p=0; p<nblocks; p++) {
             if(myid==blocks[num][p]) { myloc=p; break;}
         }

         for(int p=0; p<nblocks; p++) {
             hier::Box<NDIM> fill_box;
             int myleft = myloc-1;
             if(myleft==-1) myleft=nblocks-1;
             int myright= myloc+1;
             if(myright==nblocks) myright=0;
             if(p==myloc || p==myleft || p==myright) {
                fill_box = patches[num][p]->getBox()*
                                 center_var_scratch->getGhostBox();
                center_var_scratch->fillAll(0.0,fill_box);
             } else {
                fill_box = fill_boxes(num);
             }

             tbox::pout
                   << "\t    sing_patches = " << patches[num][p]->getBox()
                   << " fbox=" << fill_box
                   << " from block_" << blocks[num][p] << endl;
             d_math_ops->add(center_var_scratch,
                             center_var_scratch,
                             patches[num][p]
                                  ->getPatchData(d_center_var_scratch_id),
                             fill_box);
         }

      }  // end for locations of singularities.

   }  // end for types of singularities.

   for(int btype=2; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      tbox::Array< tbox::Array< tbox::Pointer <hier::Patch<NDIM> > > > patches;
      tbox::Array< tbox::Array< int > > blocks;
      tbox::Array< int > conditions;
      pgeom->getBlockPhysicalSingularBoundaryFillBoxes(fill_boxes,
                                               loc_indexes,
                                               patches,
                                               blocks,
                                               conditions,
                                               btype,
                                               patch_box, 
                                               cghost);

      const int myid = pgeom->getBlockNumber();

      for (int num=0; num < loc_indexes.size(); num++) {

         tbox::pout << "\tblock_physical_fill_box="<< fill_boxes(num)
                    << ", type=" << btype
                    << ", loc=" << loc_indexes[num]
                    << ", condition=" << conditions[num]
                    << endl;

         center_var_scratch->fillAll(0.0,fill_boxes(num));

         int myloc = -1;
         for(int p=0; p<blocks[num].size(); p++) {
             if(myid==blocks[num][p]) { myloc=p; break;}
         }

         for(int p=0; p<patches[num].size(); p++) {
             hier::Box<NDIM> fill_box;
             if(p==myloc || p==myloc-1 || p==myloc+1) {
                fill_box = patches[num][p]->getBox()*
                                 center_var_scratch->getGhostBox();
                center_var_scratch->fillAll(0.0,fill_box);
             } else {
                fill_box = fill_boxes(num);
             }

             tbox::pout 
                   << "\t    sing_patches = " << patches[num][p]->getBox()
                   << " fbox=" << fill_box
                   << " from block_" << blocks[num][p] << endl;
             d_math_ops->add(center_var_scratch,
                             center_var_scratch,
                             patches[num][p]
                                  ->getPatchData(d_center_var_scratch_id),
                             fill_box);
         }
         
      }  // end for locations of singularities.

   }  // end for types of singularities.

   const hier::Index<NDIM> ifirst=patch_box.lower();
   const hier::Index<NDIM> ilast =patch_box.upper();
   summing_(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            cghost(0),cghost(1),
            center_var_scratch->getPointer(),
            center_var_new->getPointer());

}

/*
 ********************************************************************
 * 设置物理边界条件.                                                *
 ********************************************************************
 */ 
void FDMBSingDemo::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{

   NULL_USE(fill_time);

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="SUMMING");
   #endif

   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch =
                               patch.getPatchData(d_center_var_scratch_id);

   hier::IntVector<NDIM> ghost_cells = center_var_scratch->getGhostCellWidth();

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_cells == d_oneghosts);
   TBOX_ASSERT(!center_var_scratch.isNull());
   #endif

   const tbox::Pointer< hier::FederalBlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!pgeom.isNull());
   TBOX_ASSERT(pgeom->isTouchesPhysicalBoundary());
   #endif

   const hier::Box<NDIM>& patch_box = patch.getBox();

   tbox::pout << "\nintc=" << intc_name <<",setPhysicalBoundaryConditions:" << endl;
   tbox::pout<< "block_" << pgeom->getBlockNumber()
             << ", patch_box = " << patch_box << endl;

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      tbox::Array<int> bdry_cond;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, bdry_cond,
                                          btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);
 
         tbox::pout << "\tfill_box="<< fill_box << ", type=" << btype
                    << ", loc=" << loc_indexes[num]
                    << ", condition=" << bdry_cond[num]
                    << endl;

         center_var_scratch->fillAll(0.0,fill_box);

      }

   }

}

/*
*************************************************************************
* 将FDMBSingDemo类对象的所有数据成员输出到指定的输出流.                        *
*************************************************************************
*/
void FDMBSingDemo::printClassData(ostream &os) const 
{

   os << "\nFDMBSingDemo::printClassData..." << endl;
   os << "FDMBSingDemo: this = " << (FDMBSingDemo*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::MultiblockDeformingGridGeometry<NDIM>*)d_grid_geometry << endl;

   os << "Parameters for physical problem ..." << endl;

}

/*
*************************************************************************
* 从输入数据库读取数据. 如果从输入数据库和重启动数据库中读取相同数据成  *
* 员的值.那么用输入数据库中的值覆盖重启动数据库的值.                    *
*************************************************************************
*/

void FDMBSingDemo::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
}

/*
*************************************************************************
* 将FDMBSingDemo类的数据成员写入到重启动数据库.                                *
*************************************************************************
*/

void FDMBSingDemo::putToDatabase(tbox::Pointer<tbox::Database> db)
{
}

/*
*************************************************************************
*    从重启动数据库中读取数据.                                          *
*************************************************************************
*/
void FDMBSingDemo::getFromRestart()
{
}

#endif

