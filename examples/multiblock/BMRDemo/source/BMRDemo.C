//
// 文件名:  BMRDemo.C
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 170 $
// 修改  :  $Date: 2007-06-27 08:44:22  $
// 描述  :  BMR演示: 单个网格片上的数值计算子程序.
//

#ifndef included_BMRDemo_C
#define included_BMRDemo_C

#include "BMRDemo.h"

using namespace std;

#include "BlockPatchGeometry.h"

#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/IEEE.h"
#include "tbox/Array.h"

//声明Fortran数值子程序的头文件.
#include "BMRDemoFort.h"

//解向量的分量个数
#define NEQU  (NDIM + 2) //NDIM个速度分量+ 压力 + 密度

/*
*************************************************************************
*                                                                       *
*  BMRDemo类的构造函数, 执行以下操作:                                     *
*  (1) 初始化数据成员                                                   *
*  (2) 创建定义物理模型控制方程中的数值解的变量.                        *
*  (3) 如果进行断点续算, 那么调用getFromRestart().                      *
*  (4) 调用getFromInput(), 从指定数据库中读取数据.                      *
*      该操作可能覆盖从重启动文件中读入的数据.                          *
*                                                                       *
*************************************************************************
*/
BMRDemo::BMRDemo(const string& object_name,
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
   d_grid_tool     = grid_tool;

   d_times_tag = 0;
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

/*
*************************************************************************
*    BMRDemo 类的空析构函数.	 	                                *
*************************************************************************
*/

BMRDemo::~BMRDemo() 
{ 

tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);

if(d_math_ops) delete d_math_ops;

} 

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/

void BMRDemo::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量: 坐标，中心量，结点量.
   d_coords     = new pdat::NodeVariable<NDIM,double>("coordinates", NDIM);
   d_source = new pdat::CellVariable<NDIM,double>("source", 1);
   d_center_var = new pdat::CellVariable<NDIM,double>("center_var", 1);
   d_node_var   = new pdat::NodeVariable<NDIM,double>("node_var", 1);

   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");
   d_bmr_interp     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("BMR_INTERP");

   // 存储当前时刻的变量值: (coords,current), 影像区宽度为0.
   d_coords_current_id =
         variable_db->registerVariableAndContext(d_coords,
                                                 d_current,
                                                 d_zeroghosts);

   // 注册1个影像区宽度等于2的物理量：验证几何奇异网格片的宽度.
   d_coords_scratch_id =
         variable_db->registerVariableAndContext(d_coords,
                                                 d_scratch,
                                                 d_twoghosts);

   // 存储变量的演算值: (center_var,scratch), 影像区宽度为等于1.
   d_source_id = 
         variable_db->registerVariableAndContext(d_source,
                                                 d_current,
                                                 d_zeroghosts);

   // 存储当前时刻的变量值: (center_var,current), 影像区宽度为0.
   d_center_var_current_id = 
         variable_db->registerVariableAndContext(d_center_var,
                                                 d_current,
                                                 d_zeroghosts);

   // 存储新时刻的变量值: (center_var,new), 影像区宽度为0.
   d_center_var_new_id = 
         variable_db->registerVariableAndContext(d_center_var,
                                                 d_new,
                                                 d_zeroghosts);

   // 存储变量的演算值: (center_var,new), 用于BMR插值, 影像区宽度取决于插值算子宽度.
   d_center_var_bmr_id = 
         variable_db->registerVariableAndContext(d_center_var,
                                                 d_bmr_interp,
                                                 d_oneghosts);

   // 存储变量的演算值: (center_var,scratch), 影像区宽度为等于1.
   d_center_var_scratch_id = 
         variable_db->registerVariableAndContext(d_center_var,
                                                 d_scratch,
                                                 d_oneghosts);

   // 存储当前时刻的变量值: (node_var,current), 影像区宽度为0.
   d_node_var_current_id = 
         variable_db->registerVariableAndContext(d_node_var,
                                                 d_current,
                                                 d_zeroghosts);

}

/*
*************************************************************************
*                                                                       *
* 注册JaVis数据输出器, 以采用JaVis工具对绘图文件进行后处理.           *
*                                                                       *
*************************************************************************
*/
void BMRDemo::registerJaVisDataWriter(
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

#if 1 
      d_javis_writer->registerPlotQuantity(
                              "source",
                              "SCALAR",
                              d_source_id,
                              0,1);

      d_javis_writer->registerPlotQuantity(
                              "center_var",
                              "SCALAR",
                              d_center_var_current_id,
                              0,1);

      d_javis_writer->registerPlotQuantity(
                              "node_var",
                              "SCALAR",
                              d_node_var_current_id,
                              0,1);
#endif

   }
   #endif

}
	 
/********************************************************************
*  初始化积分构件.
*********************************************************************/
void BMRDemo::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   tbox::plog << "intc_name in BMRDemo::initializeComponent= " << intc_name << endl;

   if(intc_name=="INIT") {  // 初值构件 : 为网格坐标和数值解赋初值.
        intc->registerInitPatchData(d_coords_current_id);
        intc->registerInitPatchData(d_source_id);
        intc->registerInitPatchData(d_center_var_current_id);
        intc->registerInitPatchData(d_node_var_current_id);
        intc->registerInitPatchData(d_center_var_bmr_id);
   }else if(intc_name=="TIME_STEP_SIZE") { // 步长构件：求步长.
   }else if(intc_name=="COMPUT_TEMP") { //数值构件: 求中心量的环绕累加和.
        intc->registerRefinePatchData(d_center_var_scratch_id,
                                      d_center_var_current_id);
   }else if(intc_name=="ACCEPT_SOLUTION") { //复制构件: 接收数值解.
        intc->registerCopyPatchData(d_center_var_current_id,
                                    d_center_var_new_id);
        intc->registerCopyPatchData(d_center_var_bmr_id,
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
void BMRDemo::initializePatchData( hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const bool    initial_time,
                                 const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="INIT");
   #endif

   (void) time;

   const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();

   hier::Box<NDIM> patch_box = patch.getBox();

   const hier::Index<NDIM> ifirst=patch_box.lower();
   const hier::Index<NDIM> ilast =patch_box.upper();
 
   int block_id = pgeom->getBlockNumber();

   if (initial_time) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > center_var_current = 
                               patch.getPatchData(d_center_var_current_id);
      tbox::Pointer< pdat::NodeData<NDIM,double> > node_var_current = 
                               patch.getPatchData(d_node_var_current_id);

      tbox::Pointer< pdat::CellData<NDIM,double> > source = 
                               patch.getPatchData(d_source_id);

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!coords_current.isNull());
      TBOX_ASSERT(!center_var_current.isNull());
      TBOX_ASSERT(!node_var_current.isNull());
      TBOX_ASSERT(!source.isNull());
#endif

      hier::IntVector<NDIM> ghost_cells = center_var_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ghost_cells == d_zeroghosts);
#endif

      // 生成初始网格.
      if(!d_grid_tool.isNull()) {
         d_grid_tool->generateDeformingMeshForDomain(patch,
                                                     d_coords_current_id);
      }

      // 赋初值.
      center_var_current->fillAll(0.0);
      node_var_current->fillAll(0.0);
      source->fillAll(0.0);

      if (block_id == 0 && ifirst(1) == 0) {
         init_(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               source->getPointer());
      }

   } else { // regrid time.
      tbox::plog << "BMRDemo::initializePatchData: at regrid time." << endl; 
   }

}

/*
*************************************************************************
* 初值构件: 网格结点坐标和数值解.
*************************************************************************
*/
void BMRDemo::BMRInterpolationPatchData( hier::Patch<NDIM>& patch_dst,
                                            hier::Patch<NDIM>& patch_src,
                                            const hier::IntVector<NDIM>  ratio_dst_to_src,
                                            const double  time,
                                            const bool    initial_time,
                                            const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="INIT");
   TBOX_ASSERT(initial_time==false);
   #endif

   (void) time;

   hier::Box<NDIM> patch_box_src = patch_src.getBox();
   const hier::Index<NDIM> ifirst_src=patch_box_src.lower();
   const hier::Index<NDIM> ilast_src =patch_box_src.upper();

   hier::Box<NDIM> patch_box_dst = patch_dst.getBox();
   const hier::Index<NDIM> ifirst_dst=patch_box_dst.lower();
   const hier::Index<NDIM> ilast_dst =patch_box_dst.upper();

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_bmr = 
                              patch_src.getPatchData(d_coords_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_bmr = 
                              patch_src.getPatchData(d_center_var_bmr_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > node_var_bmr = 
                              patch_src.getPatchData(d_node_var_current_id);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                              patch_dst.getPatchData(d_coords_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_current = 
                              patch_dst.getPatchData(d_center_var_current_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > node_var_current = 
                              patch_dst.getPatchData(d_node_var_current_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > source = 
                              patch_dst.getPatchData(d_source_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > source_bmr = 
                              patch_src.getPatchData(d_source_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coords_current.isNull());
   TBOX_ASSERT(!center_var_current.isNull());
   TBOX_ASSERT(!node_var_current.isNull());
   TBOX_ASSERT(!source.isNull());
   TBOX_ASSERT(!source_bmr.isNull());
   TBOX_ASSERT(!coords_bmr.isNull());
   TBOX_ASSERT(!center_var_bmr.isNull());
   TBOX_ASSERT(!node_var_bmr.isNull());
#endif
   hier::IntVector<NDIM> ghost_bmr;


/* for coordination. */
   ghost_bmr = coords_bmr->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_bmr == d_zeroghosts);
#endif

   interpolating_coord_(ifirst_src(0), ilast_src(0),
                    ifirst_src(1), ilast_src(1),
                    ghost_bmr(0), ghost_bmr(1),
                    ifirst_dst(0), ilast_dst(0),
                    ifirst_dst(1), ilast_dst(1),
                    ratio_dst_to_src(0), ratio_dst_to_src(1),
                    coords_bmr->getPointer(),
                    coords_current->getPointer());


/* for node data. */
   ghost_bmr = node_var_bmr->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_bmr == d_zeroghosts);
#endif
   interpolating_node_(ifirst_src(0), ilast_src(0),
                    ifirst_src(1), ilast_src(1),
                    ghost_bmr(0), ghost_bmr(1),
                    ifirst_dst(0), ilast_dst(0),
                    ifirst_dst(1), ilast_dst(1),
                    ratio_dst_to_src(0), ratio_dst_to_src(1),
                    node_var_bmr->getPointer(),
                    node_var_current->getPointer());

/* for cell data. */
   ghost_bmr = center_var_bmr->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_bmr == d_oneghosts);
#endif

  // fill physical boundary if needed.
   hier::IntVector<NDIM> ghost_cells;
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch;
   const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > pgeom =
                                            patch_src.getPatchGeometry();
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!pgeom.isNull());
   #endif

   bool is_touch_bdry = pgeom->isTouchesPhysicalBoundary();
   const hier::Box<NDIM>& bounding = center_var_bmr->getGhostBox();

   if (is_touch_bdry) {
   const tbox::Array<tbox::Pointer<hier::BlockPatchBoundary<NDIM> > > patch_bdry
                                       =  pgeom->getSpecialBoundaries(bounding);

   for (int pb=0; pb<patch_bdry.size(); pb++) {
       tbox::Array<tbox::Pointer<hier::Patch<NDIM> > > ghost_patches;

       if (patch_bdry[pb]->isRegularPhysicalBoundary()) {
          //tbox::pout << "patch_bdry.isRegularPhysicalBoundary(): " << endl;
          center_var_scratch = patch_src.getPatchData(d_center_var_bmr_id);
          ghost_cells = center_var_scratch->getGhostCellWidth();

          #ifdef DEBUG_CHECK_ASSERTIONS
          TBOX_ASSERT(ghost_cells == d_oneghosts);
          TBOX_ASSERT(!center_var_scratch.isNull());
          #endif
          center_var_scratch->fillAll(0.0,patch_bdry[pb]->getGhostBox(bounding));
       } else if (patch_bdry[pb]->isBlockPhysicalBoundary()) {
          //tbox::pout << "patch_bdry.isBlockPhysicalBoundary(): " << endl;
          ghost_patches = patch_bdry[pb]->getGhostPatches();
          for (int gp=0; gp<ghost_patches.size(); gp++) { 
              if (ghost_patches[gp]) {
                  center_var_scratch = ghost_patches[gp]->getPatchData(d_center_var_bmr_id);
                  #ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(ghost_cells == d_oneghosts);
                  TBOX_ASSERT(!center_var_scratch.isNull());
                  #endif
                  center_var_scratch->fillAll(0.0);
              }
          }
       }

   }
   }

   interpolating_center_(ifirst_src(0), ilast_src(0),
                    ifirst_src(1), ilast_src(1),
                    ghost_bmr(0), ghost_bmr(1),
                    ifirst_dst(0), ilast_dst(0),
                    ifirst_dst(1), ilast_dst(1),
                    ratio_dst_to_src(0), ratio_dst_to_src(1),
                    source_bmr->getPointer(),
                    source->getPointer(),
                    center_var_bmr->getPointer(),
                    center_var_current->getPointer());

    const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > pgeom_dst =
                                            patch_dst.getPatchGeometry();
    int block_id = pgeom_dst->getBlockNumber();

    source->fillAll(0.0);
    
    if (block_id == 0 && ifirst_dst(1) == 0) {
       init_(ifirst_dst(0), ilast_dst(0),
             ifirst_dst(1), ilast_dst(1),
             source->getPointer());
    }

}
/*
*************************************************************************
* 步长构件: 计算并返回网格片上的稳定时间步长.  
*************************************************************************
*/
double BMRDemo::getPatchDt(hier::Patch<NDIM>& patch,
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

   return(1.0e-1);

}

/********************************************************************
*  实现数值构件.
*********************************************************************/
void BMRDemo::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="COMPUT_TEMP");
   #endif

   if(intc_name=="COMPUT_TEMP") {
      //tbox::pout << "\nintc=" << intc_name <<",computeOnPatch:" << endl;
      advancing_temperature(patch,time,dt,initial_time);
   }

}

// 求密度的和.
void BMRDemo::advancing_temperature(hier::Patch<NDIM>& patch,
                                     const double  time,
                                     const double  dt,
                                     const bool    initial_time)
{

   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch_ghost;
   tbox::Pointer< pdat::CellData<NDIM,double> > source =
                               patch.getPatchData(d_source_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch =
                               patch.getPatchData(d_center_var_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_new =
                               patch.getPatchData(d_center_var_new_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > node_var_current =
                               patch.getPatchData(d_node_var_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!source.isNull());
   TBOX_ASSERT(!center_var_scratch.isNull());
   TBOX_ASSERT(!center_var_new.isNull());
   TBOX_ASSERT(!node_var_current.isNull());
#endif

   hier::IntVector<NDIM> cghost = center_var_scratch->getGhostCellWidth();
   hier::IntVector<NDIM> nghost = node_var_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(cghost == d_oneghosts);
   TBOX_ASSERT(nghost == d_zeroghosts);
#endif

   const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > ghost_pgeom;
   const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();

   hier::Box<NDIM> patch_box = patch.getBox();

   const hier::Index<NDIM> ifirst=patch_box.lower();
   const hier::Index<NDIM> ilast =patch_box.upper();
 
   int block_id = pgeom->getBlockNumber();
   int ghost_block_id;

   const hier::Box<NDIM>& bounding = center_var_scratch->getGhostBox();
   const tbox::Array<tbox::Pointer<hier::BlockPatchBoundary<NDIM> > > patch_bdry =
                                         pgeom->getSpecialBoundaries(bounding);


   for (int pb=0; pb<patch_bdry.size(); pb++) {
       tbox::Array<tbox::Pointer<hier::Patch<NDIM> > > ghost_patches;

       hier::Box<NDIM> fill_box;
       fill_box = patch_bdry[pb]->getGhostBox(bounding);

       ghost_patches = patch_bdry[pb]->getGhostPatches();

       int ghost_patches_size = ghost_patches.size();

       if (patch_bdry[pb]->isNode()) {

          if (!patch_bdry[pb]->isStructure()) {

             center_var_scratch->fillAll(0.0,fill_box);

             int ghost_patches_size_nonempty = 0; 
             // for each ghost patches, add to sigular fill box.
             for (int gp=0; gp<ghost_patches_size; gp++) {
                 if (ghost_patches[gp]) {
                  
                    hier::IntVector<NDIM> ratio;
                    ratio = patch_bdry[pb]->getRatio(gp);

                    ghost_block_id = patch_bdry[pb]->getBlockNumber(gp);

                    center_var_scratch_ghost = ghost_patches[gp]->getPatchData(d_center_var_scratch_id);
                    hier::Box<NDIM> ghost_box = ghost_patches[gp]->getBox();

                    hier::IntVector<NDIM> gcw_ghosts = center_var_scratch_ghost->getGhostCellWidth();

                    if (block_id != ghost_block_id) {
                        ghost_patches_size_nonempty ++;
                        const hier::Index<NDIM> ifirst_ghost=ghost_box.lower();
                        const hier::Index<NDIM> ilast_ghost =ghost_box.upper();
                        const hier::Index<NDIM> ifirst_fill=fill_box.lower();
                        const hier::Index<NDIM> ilast_fill =fill_box.upper();
                        adding_singular_(ifirst_ghost(0), ilast_ghost(0),
                                         ifirst_ghost(1), ilast_ghost(1),
                                         gcw_ghosts(0), gcw_ghosts(1),
                                         ifirst(0), ilast(0),
                                         ifirst(1), ilast(1),
                                         cghost(0),cghost(1),
                                         ifirst_fill(0), ilast_fill(0),
                                         ifirst_fill(1), ilast_fill(1),
                                         ratio(0),ratio(1),
                                         center_var_scratch_ghost->getPointer(),
                                         center_var_scratch->getPointer());
 
                    }

                 }
             } // end for ghost patches

             if (ghost_patches_size_nonempty != 0) { 
                double scaling = 1.0/(ghost_patches_size_nonempty);
                //tbox::pout << " scaling = " << scaling << endl;
                d_math_ops->scale(center_var_scratch,
                                  scaling,
                                  center_var_scratch,
                                  fill_box);
             }

          } else if (patch_bdry[pb]->isOne2Many()){ // 一对多情形, 是否该有个特定的判断. 

            //tbox::pout << "patch_bdry.isStructure and One2Many: " << endl;
            
            // in this case, fill_box may be too large, special treatment performed in fortran subroutine. 
            const hier::Index<NDIM> ifirst_fill_0=fill_box.lower();
            const hier::Index<NDIM> ilast_fill_0 =fill_box.upper();
            int size_0 = ilast_fill_0(0)-ifirst_fill_0(0)+1;
            int size_1 = ilast_fill_0(1)-ifirst_fill_0(1)+1;
            hier::IntVector<NDIM> size_reduce;

            size_reduce(0) = cghost(0) - size_0;
            size_reduce(1) = cghost(1) - size_1;

            fill_box.growUpper(size_reduce);

            center_var_scratch->fillAll(0.0,fill_box);
           
            int ghost_patches_size_nonempty = 0; 
            for (int gp=0; gp<ghost_patches.size(); gp++) {
                if (ghost_patches[gp]) {
                
                    ghost_patches_size_nonempty++;

                    hier::IntVector<NDIM> ratio;
                    ratio = patch_bdry[pb]->getRatio(gp);

                    ghost_block_id = patch_bdry[pb]->getBlockNumber(gp);

                    center_var_scratch_ghost = ghost_patches[gp]->getPatchData(d_center_var_scratch_id);
                    hier::Box<NDIM> ghost_box = ghost_patches[gp]->getBox();

                    hier::IntVector<NDIM> gcw_ghosts = center_var_scratch_ghost->getGhostCellWidth();

                    const hier::Index<NDIM> ifirst_ghost=ghost_box.lower();
                    const hier::Index<NDIM> ilast_ghost =ghost_box.upper();
                    const hier::Index<NDIM> ifirst_fill=fill_box.lower();
                    const hier::Index<NDIM> ilast_fill =fill_box.upper();
                 
                    filling_one2many_node_(ifirst_ghost(0), ilast_ghost(0),
                                     ifirst_ghost(1), ilast_ghost(1),
                                     gcw_ghosts(0), gcw_ghosts(1),
                                     ifirst(0), ilast(0),
                                     ifirst(1), ilast(1),
                                     cghost(0),cghost(1),
                                     ifirst_fill(0), ilast_fill(0),
                                     ifirst_fill(1), ilast_fill(1),
                                     ratio(0),ratio(1),
                                     center_var_scratch_ghost->getPointer(),
                                     center_var_scratch->getPointer());

                }
            }

            if (ghost_patches_size_nonempty != 0) {
               double scaling = 1.0/(ghost_patches_size_nonempty);
               d_math_ops->scale(center_var_scratch,
                                 scaling,
                                 center_var_scratch,
                                 fill_box);
            }

         } else {
           //tbox::plog << "patch_bdry.isStructure and One2One: no need treatment." << endl;
         } 
       }

       if (patch_bdry[pb]->isEdge()) {
          if (patch_bdry[pb]->isOne2Many()) {
             for (int gp=0; gp<ghost_patches.size(); gp++) {
                 if (ghost_patches[gp]) {

                    hier::IntVector<NDIM> ratio;
                    ratio = patch_bdry[pb]->getRatio(gp);

                    ghost_block_id = patch_bdry[pb]->getBlockNumber(gp);

                    center_var_scratch_ghost = ghost_patches[gp]->getPatchData(d_center_var_scratch_id);
                    hier::Box<NDIM> ghost_box = ghost_patches[gp]->getBox();

                    hier::IntVector<NDIM> gcw_ghosts = center_var_scratch_ghost->getGhostCellWidth();

                    const hier::Index<NDIM> ifirst_ghost=ghost_box.lower();
                    const hier::Index<NDIM> ilast_ghost =ghost_box.upper();
                    const hier::Index<NDIM> ifirst_fill=fill_box.lower();
                    const hier::Index<NDIM> ilast_fill =fill_box.upper();
                    filling_one2many_(ifirst_ghost(0), ilast_ghost(0),
                                     ifirst_ghost(1), ilast_ghost(1),
                                     gcw_ghosts(0), gcw_ghosts(1),
                                     ifirst(0), ilast(0),
                                     ifirst(1), ilast(1),
                                     cghost(0),cghost(1),
                                     ifirst_fill(0), ilast_fill(0),
                                     ifirst_fill(1), ilast_fill(1),
                                     ratio(0),ratio(1),
                                     center_var_scratch_ghost->getPointer(),
                                     center_var_scratch->getPointer());

                 }
             }
          } else { // 在物理边界中处理.
             //tbox::plog << "patch_bdry.isOne2One: treat in PhysicalBoundary." << endl;
          }
       }

   }

   advancing_(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            cghost(0),cghost(1),
            dt,
            node_var_current->getPointer(),
            source->getPointer(),
            center_var_scratch->getPointer(),
            center_var_new->getPointer());

}

/*
 ********************************************************************
 * 设置物理边界条件.                                                *
 ********************************************************************
 */ 
void BMRDemo::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{

   NULL_USE(fill_time);

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(intc_name=="COMPUT_TEMP");
   #endif

   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_current =
                               patch.getPatchData(d_center_var_current_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > center_var_scratch;
   hier::IntVector<NDIM> ghost_cells;

#if 1 
   center_var_scratch = patch.getPatchData(d_center_var_scratch_id);
   ghost_cells = center_var_scratch->getGhostCellWidth();

   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_cells == d_oneghosts);
   TBOX_ASSERT(!center_var_current.isNull());
   TBOX_ASSERT(!center_var_scratch.isNull());
   #endif
#endif

   const tbox::Pointer< hier::AdvancedBlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!pgeom.isNull());
   TBOX_ASSERT(pgeom->isTouchesPhysicalBoundary());
   #endif

   const hier::Box<NDIM>& patch_box = patch.getBox();

   const hier::Index<NDIM> ifirst=patch_box.lower();
   const hier::Index<NDIM> ilast =patch_box.upper();
 
   int block_id = pgeom->getBlockNumber();

   const hier::Box<NDIM>& bounding = center_var_scratch->getGhostBox();

   const tbox::Array<tbox::Pointer<hier::BlockPatchBoundary<NDIM> > > patch_bdry
                                       =  pgeom->getSpecialBoundaries(bounding);

   for (int pb=0; pb<patch_bdry.size(); pb++) {
       tbox::Array<tbox::Pointer<hier::Patch<NDIM> > > ghost_patches;

       hier::Box<NDIM> fill_box;
       fill_box = patch_bdry[pb]->getGhostBox(bounding);

       if (block_id == 0 || block_id == 1 || block_id == 2 || block_id == 3 || block_id == 4 || block_id == 5) {

       const hier::Index<NDIM> ifirst_f=fill_box.lower();
       const hier::Index<NDIM> ilast_f =fill_box.upper();

       filling_phybdr_(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            ghost_cells(0),ghost_cells(1),
            ifirst_f(0), ilast_f(0),
            ifirst_f(1), ilast_f(1),
            center_var_current->getPointer(),
            center_var_scratch->getPointer());

       } else { 

       if (patch_bdry[pb]->isRegularPhysicalBoundary()) {
          //tbox::pout << "patch_bdry.isRegularPhysicalBoundary(): " << endl;
          center_var_scratch = patch.getPatchData(d_center_var_scratch_id);
          ghost_cells = center_var_scratch->getGhostCellWidth();

          #ifdef DEBUG_CHECK_ASSERTIONS
          TBOX_ASSERT(ghost_cells == d_oneghosts);
          TBOX_ASSERT(!center_var_scratch.isNull());
          #endif
          center_var_scratch->fillAll(0.0,fill_box);
       } else if (patch_bdry[pb]->isBlockPhysicalBoundary()) {
          ghost_patches = patch_bdry[pb]->getGhostPatches();
          for (int gp=0; gp<ghost_patches.size(); gp++) { 
              if (ghost_patches[gp]) {
                  center_var_scratch = ghost_patches[gp]->getPatchData(d_center_var_scratch_id);
                  #ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(ghost_cells == d_oneghosts);
                  TBOX_ASSERT(!center_var_scratch.isNull());
                  #endif
                  center_var_scratch->fillAll(0.0);
              }
          }
       }

       }
   }

}

/*
*************************************************************************
* 将BMRDemo类对象的所有数据成员输出到指定的输出流.                        *
*************************************************************************
*/
void BMRDemo::printClassData(ostream &os) const 
{

   os << "\nBMRDemo::printClassData..." << endl;
   os << "BMRDemo: this = " << (BMRDemo*)this << endl;
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

void BMRDemo::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
}

/*
*************************************************************************
* 将BMRDemo类的数据成员写入到重启动数据库.                                *
*************************************************************************
*/

void BMRDemo::putToDatabase(tbox::Pointer<tbox::Database> db)
{
}

/*
*************************************************************************
*    从重启动数据库中读取数据.                                          *
*************************************************************************
*/
void BMRDemo::getFromRestart()
{
}

#endif

