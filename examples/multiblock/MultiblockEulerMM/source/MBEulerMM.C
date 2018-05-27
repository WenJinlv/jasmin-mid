//
// 文件名:  MBEulerMM.C
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 170 $
// 修改  :  $Date: 2007-06-27 08:44:22  $
// 描述  :  Euler方程组：移动网格方法在单网格片上的子程序.
//

#ifndef included_MBEulerMM_C
#define included_MBEulerMM_C

#include "MBEulerMM.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

using namespace std;

#include "BlockPatchGeometry.h"

#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/IEEE.h"

//声明Fortran数值子程序的头文件.
#include "MBEulerMMFort.h"

//解向量的分量个数
#define NEQU  (NDIM + 2) //NDIM个速度分量+ 压力 + 密度

/*
*************************************************************************
*                                                                       *
*  MBEulerMM类的构造函数, 执行以下操作:                                     *
*  (1) 初始化数据成员                                                   *
*  (2) 创建定义物理模型控制方程中的数值解的变量.                        *
*  (3) 如果进行断点续算, 那么调用getFromRestart().                      *
*  (4) 调用getFromInput(), 从指定数据库中读取数据.                      *
*      该操作可能覆盖从重启动文件中读入的数据.                          *
*                                                                       *
*************************************************************************
*/
MBEulerMM::MBEulerMM(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
         tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
   assert(!grid_tool.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometry = grid_geom;
   d_grid_tool = grid_tool;

   // 缺省值.
   d_gamma = 1.4;
   d_alpha = 0.0;
   d_beta  = 0.9;

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

   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();

}

/*
*************************************************************************
*    MBEulerMM 类的空析构函数.	 	                                *
*************************************************************************
*/

MBEulerMM::~MBEulerMM() 
{ 

tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
} 

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/

void MBEulerMM::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量. NEQU: 0表示密度, 1,2表示动量, 3表示能量密度.
   d_coords   = new pdat::NodeVariable<NDIM,double>("coordinates", NDIM);
   d_solution = new pdat::CellVariable<NDIM,double>("solution", NEQU);
   d_flux     = new pdat::SideVariable<NDIM,double>("flux", NEQU);
   d_wc       = new pdat::CellVariable<NDIM,double>("wc",1);
   d_dux      = new pdat::CellVariable<NDIM,double>("dux",4);
   d_duy      = new pdat::CellVariable<NDIM,double>("duy",4);

   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");
   d_move = hier::VariableDatabase<NDIM>::getDatabase()->getContext("MOVE");

   // 存储当前时刻的变量值: (coords,current), 影像区宽度为0.
   d_coords_current_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_current,
                                                 d_zeroghosts);

   // 存储新时刻的变量值: (coords,new), 影像区宽度为0.
   d_coords_new_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_new,
                                                 d_zeroghosts);

   // 存储变量的演算值: (coords,scratch), 影像区宽度等于格式宽度.
   d_coords_scratch_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_scratch,
                                                 d_oneghosts);
   // 存储变量的演算值: (coords,move), 影像区宽度等于格式宽度.
   d_coords_move_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_move,
                                                 d_oneghosts);

   // 存储当前时刻的变量值: (solution,current), 影像区宽度为0.
   d_solution_current_id = 
         variable_db->registerVariableAndContext(d_solution,
                                                 d_current,
                                                 d_zeroghosts);

   // 存储新时刻的变量值: (solution,new), 影像区宽度为0.
   d_solution_new_id = 
         variable_db->registerVariableAndContext(d_solution,
                                                 d_new,
                                                 d_zeroghosts);

   // 存储变量的演算值: (solution,scratch), 影像区宽度为等于1.
   d_solution_scratch_id = 
         variable_db->registerVariableAndContext(d_solution,
                                                 d_scratch,
                                                 d_oneghosts);

   // 存储新时刻的通量值: (flux,new), 影像区宽度等于0.
   d_flux_new_id = 
          variable_db->registerVariableAndContext(d_flux,
                                                  d_new,
                                                  d_zeroghosts);

   // 存储新时刻的变量值: (wc,new), 影像区宽度为0.
   d_wc_new_id = 
         variable_db->registerVariableAndContext(d_wc,
                                                 d_new,
                                                 d_zeroghosts);

   // 存储变量的演算值: (wc,scratch), 影像区宽度为等于1.
   d_wc_scratch_id = 
         variable_db->registerVariableAndContext(d_wc,
                                                 d_scratch,
                                                 d_oneghosts);

   // 存储变量的演算值: (dux,scratch), 影像区宽度为等于1.
   d_dux_scratch_id = 
         variable_db->registerVariableAndContext(d_dux,
                                                 d_scratch,
                                                 d_oneghosts);

   // 存储变量的演算值: (duy,scratch), 影像区宽度为等于1.
   d_duy_scratch_id = 
         variable_db->registerVariableAndContext(d_duy,
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
void MBEulerMM::registerDataWriter(
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer,
   tbox::Pointer<appu::DomainDataWriter<NDIM> > domain_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(javis_writer.isNull()));
   assert(domain_writer);
#endif
   d_javis_writer = javis_writer;

   // 注册可视化量.
   #ifdef HAVE_HDF5

   if (!(d_javis_writer.isNull())) {

      d_javis_writer->registerNodeCoordinates(d_coords_current_id);

      d_javis_writer->registerPlotQuantity(
                              "density",
                              "SCALAR",
                              d_solution_current_id,
                              0,1);

      d_javis_writer->registerPlotQuantity(
                              "momentum",
                              "VECTOR",
                              d_solution_current_id,
                              1,2);

      d_javis_writer->registerPlotQuantity(
                              "energy",
                              "SCALAR",
                              d_solution_current_id,
                              3,1);
      
   }
   #endif

   if (domain_writer){
      domain_writer->registerQuantity(d_coords_current_id);
      domain_writer->registerQuantity(d_solution_current_id);
   }

}
	 
/********************************************************************
*  初始化积分构件.
*********************************************************************/
void MBEulerMM::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件 : 为网格坐标和数值解赋初值.
        intc->registerInitPatchData(d_coords_current_id);
        intc->registerInitPatchData(d_solution_current_id);
   }else if(intc_name=="TIME_STEP_SIZE") { //步长构件: 求时间步长.

   }else if(intc_name=="COPY_PREMOVING") { //复制构件: 当前值赋值给新值.
        intc->registerCopyPatchData(d_coords_new_id,
                                d_coords_current_id);
        intc->registerCopyPatchData(d_solution_new_id,
                                d_solution_current_id);
   }else if(intc_name=="COMPUTE_WCNEW") { //数值构件: 求控制函数.
        intc->registerRefinePatchData(d_solution_scratch_id,
                                      d_solution_new_id);
        intc->registerRefinePatchData(d_coords_scratch_id,
                                      d_coords_new_id);
   }else if(   intc_name=="WCNEW_ITERATION" 
            || intc_name=="MOVING_MESH") { //数值构件: 迭代求控制函数, 移动网格.
        intc->registerRefinePatchData(d_wc_scratch_id,
                                      d_wc_new_id);
   }else if(intc_name=="MOVING_GRAD") { // 数值构件: 移动网格前, 计算梯度.

   }else if(intc_name=="MOVING_SOLUTION") { // 数值构件: 移动网格.
        intc->registerRefinePatchData(d_coords_move_id,
                                      d_coords_new_id);
        intc->registerRefinePatchData(d_dux_scratch_id,
                                      d_dux_scratch_id);
        intc->registerRefinePatchData(d_duy_scratch_id,
                                      d_duy_scratch_id);
   }else if(intc_name=="COMPUTE_GRAD") { // 数值构件: 求Euler方程, 计算梯度.
        intc->registerRefinePatchData(d_coords_scratch_id,
                                      d_coords_new_id);
        intc->registerRefinePatchData(d_solution_scratch_id,
                                      d_solution_new_id);
   }else if(intc_name=="SOLVE_EQUATION") { // 数值构件: 求Euler方程.
        intc->registerRefinePatchData(d_dux_scratch_id,
                                      d_dux_scratch_id);
        intc->registerRefinePatchData(d_duy_scratch_id,
                                      d_duy_scratch_id);
   }else if(intc_name=="ITER_ERROR") { 

   }else if(intc_name=="RESET_SOLUTION") { // 复制构件: 接收数值解.
        intc->registerCopyPatchData(d_coords_current_id,
                                    d_coords_new_id);
        intc->registerCopyPatchData(d_solution_current_id,
                                    d_solution_new_id);
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 内存构件: 为新值数据片调度内存空间.
        intc->registerPatchData(d_coords_new_id);
        intc->registerPatchData(d_solution_new_id);
        intc->registerPatchData(d_flux_new_id);    
        intc->registerPatchData(d_wc_new_id);
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 内存构件: 为演算数据片调度内存空间.
        intc->registerPatchData(d_coords_scratch_id);
        intc->registerPatchData(d_coords_move_id);
        intc->registerPatchData(d_solution_scratch_id);
        intc->registerPatchData(d_wc_scratch_id);
        intc->registerPatchData(d_dux_scratch_id);
        intc->registerPatchData(d_duy_scratch_id);
   }else if(intc_name=="REMAP"){
        intc->registerPatchData(d_coords_current_id);
        intc->registerPatchData(d_solution_current_id);
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
void MBEulerMM::initializePatchData( hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const bool    initial_time,
                                 const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT");
   #endif

   (void) time;

   if (initial_time) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > solution_current = 
                               patch.getPatchData(d_solution_current_id);


#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!coords_current.isNull());
      assert(!solution_current.isNull());
#endif

      hier::IntVector<NDIM> ghost_cells = solution_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(ghost_cells == d_zeroghosts);
#endif

      // 生成初始网格.
      if(!d_grid_tool.isNull()) {
         d_grid_tool->generateDeformingMeshForDomain(patch,
                                                     d_coords_current_id);
      }

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      // 数值解.
      tbox::Pointer< hier::BlockPatchGeometry<NDIM> > bgeom = 
                                          patch.getPatchGeometry();
      setinitialvalue_(ifirst(0),ilast(0),
                       ifirst(1),ilast(1),
                       ghost_cells(0), ghost_cells(1),
                       bgeom->getBlockNumber(),
                       0,
                       d_gamma,
                       solution_current->getPointer());

   }

}

/*
*************************************************************************
* 步长构件: 计算并返回网格片上的稳定时间步长.  
*************************************************************************
*/
double MBEulerMM::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TIME_STEP_SIZE");
#endif

   NULL_USE(flag_last_dt);
   NULL_USE(last_dt);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_current = 
                               patch.getPatchData(d_solution_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_current.isNull());
   assert(!solution_current.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = solution_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_zeroghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   double stabdt;
   stabledt_(ifirst(0), ilast(0),
             ifirst(1), ilast(1),
             ghost_cells(0), ghost_cells(1),
             d_gamma,
             coords_current->getPointer(0),
             coords_current->getPointer(1),
             solution_current->getPointer(),
             stabdt);

   return stabdt;

}

/********************************************************************
*  实现数值构件.
*********************************************************************/
void MBEulerMM::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{

   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(   intc_name=="SOLVE_EQUATION" || intc_name=="COMPUTE_GRAD"
          || intc_name=="COMPUTE_WCNEW"  || intc_name=="WCNEW_ITERATION"
          || intc_name=="MOVING_MESH"    || intc_name=="MOVING_GRAD"
          || intc_name=="MOVING_SOLUTION");
   #endif

   if(intc_name == "SOLVE_EQUATION") {
      solveEulerEquation(patch,time,dt,initial_time);
   } else if (intc_name == "COMPUTE_GRAD") {
      computeGradient(patch,time,dt,initial_time);
   } else if (intc_name == "COMPUTE_WCNEW") {
      computeWCNEW(patch,time,dt,initial_time);
   } else if (intc_name == "WCNEW_ITERATION") {
      wcnewIteration(patch,time,dt,initial_time);
   } else if (intc_name == "MOVING_MESH") {
      moveMesh(patch,time,dt,initial_time);
   } else if (intc_name == "MOVING_GRAD") {
      computeGradient(patch,time,dt,initial_time);
   } else if (intc_name == "MOVING_SOLUTION") {
      moveSolution(patch,time,dt,initial_time);
   } 

}

/********************************************************************
*  数值构件: 求Euler方程.
*********************************************************************/
void MBEulerMM::solveEulerEquation(hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const double  dt,
                                 const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch= 
                               patch.getPatchData(d_coords_scratch_id);
   tbox::Pointer< pdat::SideData<NDIM,double> > flux_new = 
                               patch.getPatchData(d_flux_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_scratch = 
                               patch.getPatchData(d_solution_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_new = 
                               patch.getPatchData(d_solution_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > dux_scratch = 
                               patch.getPatchData(d_dux_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > duy_scratch = 
                               patch.getPatchData(d_duy_scratch_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_scratch.isNull());
   assert(!solution_scratch.isNull());
   assert(!solution_new.isNull());
   assert(!flux_new.isNull());
   assert(!dux_scratch.isNull());
   assert(!duy_scratch.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = solution_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   solvepde_(ifirst(0), ilast(0),
             ifirst(1), ilast(1),
             ghost_cells(0), ghost_cells(1),
             dt,d_gamma,
             coords_scratch->getPointer(0),
             coords_scratch->getPointer(1),
             solution_scratch->getPointer(),
             dux_scratch->getPointer(),
             duy_scratch->getPointer(),
             flux_new->getPointer(0),
             flux_new->getPointer(1),
             solution_new->getPointer());

}

/********************************************************************
*  数值构件: 求梯度.
*********************************************************************/
void MBEulerMM::computeGradient(hier::Patch<NDIM>& patch,
                              const double  time,
                              const double  dt,
                              const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch= 
                               patch.getPatchData(d_coords_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_scratch = 
                               patch.getPatchData(d_solution_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > dux_scratch = 
                               patch.getPatchData(d_dux_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > duy_scratch = 
                               patch.getPatchData(d_duy_scratch_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_scratch.isNull());
   assert(!solution_scratch.isNull());
   assert(!dux_scratch.isNull());
   assert(!duy_scratch.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = solution_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   computegradient_(ifirst(0), ilast(0),
                    ifirst(1), ilast(1),
                    ghost_cells(0), ghost_cells(1),
                    coords_scratch->getPointer(0),
                    coords_scratch->getPointer(1),
                    solution_scratch->getPointer(),
                    dux_scratch->getPointer(),
                    duy_scratch->getPointer());

}

/********************************************************************
*  数值构件: 计算控制函数.
*********************************************************************/
void MBEulerMM::computeWCNEW(hier::Patch<NDIM>& patch,
                           const double  time,
                           const double  dt,
                           const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch= 
                               patch.getPatchData(d_coords_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_scratch = 
                               patch.getPatchData(d_solution_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > wc_new = 
                               patch.getPatchData(d_wc_new_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_scratch.isNull());
   assert(!solution_scratch.isNull());
   assert(!wc_new.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = solution_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   double xdx = (d_xhi[0]-d_xlo[0])/d_ncells(0);
   double ydy = (d_xhi[1]-d_xlo[1])/d_ncells(1);

   computewcnew_(ifirst(0), ilast(0),
                 ifirst(1), ilast(1),
                 ghost_cells(0), ghost_cells(1),
                 d_alpha, d_beta, xdx, ydy,
                 coords_scratch->getPointer(0),
                 coords_scratch->getPointer(1),
                 solution_scratch->getPointer(),
                 wc_new->getPointer());

}


/********************************************************************
*  数值构件: 迭代求解控制函数.
*********************************************************************/
void MBEulerMM::wcnewIteration(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > wc_new = 
                               patch.getPatchData(d_wc_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > wc_scratch = 
                               patch.getPatchData(d_wc_scratch_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!wc_new.isNull());
   assert(!wc_scratch.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = wc_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   wcnewiteration_(ifirst(0), ilast(0),
                   ifirst(1), ilast(1),
                   ghost_cells(0), ghost_cells(1),
                   wc_scratch->getPointer(),
                   wc_new->getPointer());

}

/********************************************************************
*  数值构件: 移动网格.
*********************************************************************/
void MBEulerMM::moveMesh(hier::Patch<NDIM>& patch,
                       const double  time,
                       const double  dt,
                       const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch= 
                               patch.getPatchData(d_coords_scratch_id);      
   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_new= 
                               patch.getPatchData(d_coords_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > wc_scratch = 
                               patch.getPatchData(d_wc_scratch_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_scratch.isNull());
   assert(!coords_new.isNull());
   assert(!wc_scratch.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = coords_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   movemesh_(ifirst(0), ilast(0),
             ifirst(1), ilast(1),
             ghost_cells(0), ghost_cells(1),
             coords_scratch->getPointer(0),
             coords_scratch->getPointer(1),
             wc_scratch->getPointer(),
             coords_new->getPointer(0),
             coords_new->getPointer(1));

   const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom =
                                       patch.getPatchGeometry();
   if(pgeom->isTouchesPhysicalBoundary()) {
      const hier::Box<NDIM>& patch_box = patch.getBox();

      const int btype = 1;
      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, 
                                          btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {
             hier::Box<NDIM> fill_box = fill_boxes(num);
             const hier::Index<NDIM> igfirst=fill_box.lower();
             const hier::Index<NDIM> iglast =fill_box.upper();
             postprocessbdrymesh_(ifirst(0), ilast(0),
                                  ifirst(1), ilast(1),
                                  ghost_cells(0), ghost_cells(1),
                                  igfirst(0),iglast(0),
                                  igfirst(1),iglast(1),
                                  btype, loc_indexes[num], 
                                  coords_scratch->getPointer(0),
                                  coords_scratch->getPointer(1),
                                  coords_new->getPointer(0),
                                  coords_new->getPointer(1));
      }

   }

}

/********************************************************************
*  数值构件: 移动网格之后, 修正数值解.
*********************************************************************/
void MBEulerMM::moveSolution(hier::Patch<NDIM>& patch,
                           const double  time,
                           const double  dt,
                           const bool    initial_time)
{

   NULL_USE(time);
   NULL_USE(initial_time);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch= 
                               patch.getPatchData(d_coords_scratch_id);      
   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_move = 
                               patch.getPatchData(d_coords_move_id); 
   tbox::Pointer< pdat::CellData<NDIM,double> > dux_scratch = 
                               patch.getPatchData(d_dux_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > duy_scratch = 
                               patch.getPatchData(d_duy_scratch_id);
   tbox::Pointer< pdat::SideData<NDIM,double> > flux_new = 
                               patch.getPatchData(d_flux_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_scratch = 
                               patch.getPatchData(d_solution_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_new = 
                               patch.getPatchData(d_solution_new_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_scratch.isNull());
   assert(!coords_move.isNull());
   assert(!dux_scratch.isNull());
   assert(!duy_scratch.isNull());
   assert(!flux_new.isNull());
   assert(!solution_scratch.isNull());
   assert(!solution_new.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = coords_scratch->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_oneghosts);
#endif

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   movesolution_(ifirst(0), ilast(0),
                   ifirst(1), ilast(1),
                   ghost_cells(0), ghost_cells(1),
                   coords_scratch->getPointer(0),
                   coords_scratch->getPointer(1),
                   coords_move->getPointer(0),
                   coords_move->getPointer(1),
                   solution_scratch->getPointer(),
                   solution_new->getPointer(),
                   dux_scratch->getPointer(),
                   duy_scratch->getPointer(),
                   flux_new->getPointer(0),
                   flux_new->getPointer(1));

}

/********************************************************************
 * 计算迭代误差.
 ********************************************************************/
void MBEulerMM::reduceOnPatch( 
                    double *error,
                    int     len,
                    hier::Patch<NDIM> &patch,
                    const double time,
                    const double dt,
                    const string& intc_name )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="ITER_ERROR");
#endif
              tbox::Pointer< pdat::NodeData<NDIM,double> > coords_last
                              = patch.getPatchData(d_coords_scratch_id);
              tbox::Pointer< pdat::NodeData<NDIM,double> > coords_new
                              = patch.getPatchData(d_coords_move_id);

              #ifdef DEBUG_CHECK_ASSERTIONS
              TBOX_ASSERT(coords_last->getGhostCellWidth()
                                  == hier::IntVector<NDIM>(1));       
              TBOX_ASSERT(coords_new->getGhostCellWidth() 
                                  == hier::IntVector<NDIM>(1));
              #endif

              double* clast = coords_last->getPointer();
              double* cnew  = coords_new->getPointer();

              int ncx = patch.getBox().numberCells(0)+3;
              int ncy = patch.getBox().numberCells(1)+3;

              for(int j=0; j<ncy; j++) {
                  for(int i=0; i<ncx; i++) {
                      int    id  = j*ncx+i;
                      double val = cnew[id]- clast[id];
                      if(abs(clast[id])>1.e-6) val /=  clast[id];
                      val = abs(val);
                      if(val > *error) *error = val; 
                  }
              }

} 
/*
 ********************************************************************
 * 设置物理边界条件.                                                *
 ********************************************************************
 */ 
void MBEulerMM::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{

   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(   intc_name=="SOLVE_EQUATION" || intc_name=="COMPUTE_GRAD"
          || intc_name=="COMPUTE_WCNEW"  || intc_name=="WCNEW_ITERATION"
          || intc_name=="MOVING_MESH"    || intc_name=="MOVING_GRAD"
          || intc_name=="MOVING_SOLUTION"  );
   #endif

   NULL_USE(fill_time);

   hier::IntVector<NDIM> ghost_cells(0);
   tbox::Pointer< pdat::CellData<NDIM,double> > dux ;
   tbox::Pointer< pdat::CellData<NDIM,double> > duy ; 
   tbox::Pointer< pdat::CellData<NDIM,double> > solution_scratch ; 
   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_scratch ;
   tbox::Pointer< pdat::CellData<NDIM,double> > wc_scratch ;
   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_move ;

   if(intc_name=="SOLVE_EQUATION") {

      dux = patch.getPatchData(d_dux_scratch_id);
      duy = patch.getPatchData(d_duy_scratch_id);

      ghost_cells = dux->getGhostCellWidth();

      #ifdef  DEBUG_CHECK_ASSERTIONS
      assert(!dux.isNull());
      assert(!duy.isNull());
      assert(ghost_cells == d_oneghosts);
      #endif

   } else if(intc_name=="COMPUTE_GRAD" || intc_name=="COMPUTE_WCNEW"  ) {

      solution_scratch = patch.getPatchData(d_solution_scratch_id);
      coords_scratch   = patch.getPatchData(d_coords_scratch_id);

      ghost_cells = solution_scratch->getGhostCellWidth();

      #ifdef DEBUG_CHECK_ASSERTIONS
      assert(!solution_scratch.isNull());
      assert(!coords_scratch.isNull());
      assert(ghost_cells == d_oneghosts);
      #endif

   } else if(intc_name == "WCNEW_ITERATION" || intc_name == "MOVING_MESH") {

      wc_scratch = patch.getPatchData(d_wc_scratch_id);

      ghost_cells = wc_scratch->getGhostCellWidth();
 
      #ifdef DEBUG_CHECK_ASSERTIONS
      assert(!wc_scratch.isNull());
      assert(ghost_cells == d_oneghosts);
      #endif

   } else if(intc_name == "MOVING_SOLUTION") {

      coords_move = patch.getPatchData(d_coords_move_id);
      dux = patch.getPatchData(d_dux_scratch_id);
      duy = patch.getPatchData(d_duy_scratch_id);

      ghost_cells = coords_move->getGhostCellWidth();
 
      #ifdef DEBUG_CHECK_ASSERTIONS
      assert(!coords_move.isNull());
      assert(ghost_cells == d_oneghosts);
      #endif

   }

   const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();
   #ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!pgeom.isNull());
   TBOX_ASSERT(pgeom->isTouchesPhysicalBoundary());
   #endif

   const hier::Box<NDIM>& patch_box = patch.getBox();

//   tbox::pout<< "block_" << pgeom->getBlockNumber()
//             << ", patch_box = " << patch_box << endl;

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      tbox::Array<int> bdry_cond;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes,bdry_cond,
                                          btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);

 //        tbox::pout << "\tfill_box="<< fill_box << ", type=" << btype
 //                   << ", loc=" << loc_indexes[num] 
 //                   << ", condition=" << bdry_cond[num]
 //                   << endl;

         if(intc_name=="SOLVE_EQUATION") {

            dux->fillAll(0.0,fill_box);
            duy->fillAll(0.0,fill_box);

         } else if(intc_name=="COMPUTE_GRAD"|| intc_name=="COMPUTE_WCNEW") {
            
           const hier::Index<NDIM> ifirst =patch.getBox().lower();
           const hier::Index<NDIM> ilast  =patch.getBox().upper();
           const hier::Index<NDIM> igfirst=fill_box.lower();
           const hier::Index<NDIM> iglast =fill_box.upper();

           double xdx = (d_xhi[0]-d_xlo[0])/d_ncells(0);
           double ydy = (d_xhi[1]-d_xlo[1])/d_ncells(1);

           setphysbdryfornodes_(ifirst(0), ilast(0),
                                ifirst(1), ilast(1),
                                ghost_cells(0), ghost_cells(1),
                                d_xlo[0],d_xhi[0],d_xlo[1],d_xhi[1],xdx,ydy,
                                igfirst(0), iglast(0),
                                igfirst(1), iglast(1),
                                btype, loc_indexes[num], 
                                coords_scratch->getPointer(0),
                                coords_scratch->getPointer(1));

           const int lmax = NDIM+2;
           setphysbdryforcells_(ifirst(0), ilast(0),
                                ifirst(1), ilast(1),
                                ghost_cells(0), ghost_cells(1),
                                igfirst(0), iglast(0),
                                igfirst(1), iglast(1),
                                btype, loc_indexes[num], lmax,
                                solution_scratch->getPointer());

         } else if(intc_name == "WCNEW_ITERATION" || 
                   intc_name == "MOVING_MESH") {

           const hier::Index<NDIM> ifirst =patch.getBox().lower();
           const hier::Index<NDIM> ilast  =patch.getBox().upper();
           const hier::Index<NDIM> igfirst=fill_box.lower();
           const hier::Index<NDIM> iglast =fill_box.upper();

           const int lmax = 1;
           setphysbdryforcells_(ifirst(0), ilast(0),
                                ifirst(1), ilast(1),
                                ghost_cells(0), ghost_cells(1),
                                igfirst(0), iglast(0),
                                igfirst(1), iglast(1),
                                btype, loc_indexes[num], lmax,
                                wc_scratch->getPointer());

         } else if(intc_name == "MOVING_SOLUTION") {

           const hier::Index<NDIM> ifirst =patch.getBox().lower();
           const hier::Index<NDIM> ilast  =patch.getBox().upper();
           const hier::Index<NDIM> igfirst=fill_box.lower();
           const hier::Index<NDIM> iglast =fill_box.upper();

           double xdx = (d_xhi[0]-d_xlo[0])/d_ncells(0);
           double ydy = (d_xhi[1]-d_xlo[1])/d_ncells(1);

           setphysbdryfornodes_(ifirst(0), ilast(0),
                                ifirst(1), ilast(1),
                                ghost_cells(0), ghost_cells(1),
                                d_xlo[0],d_xhi[0],d_xlo[1],d_xhi[1],xdx,ydy,
                                igfirst(0), iglast(0),
                                igfirst(1), iglast(1),
                                btype, loc_indexes[num], 
                                coords_move->getPointer(0),
                                coords_move->getPointer(1));

           dux->fillAll(0.0,fill_box);
           duy->fillAll(0.0,fill_box);
         }

      }

   }

}

/*
*************************************************************************
* 将MBEulerMM类对象的所有数据成员输出到指定的输出流.                        *
*************************************************************************
*/
void MBEulerMM::printClassData(ostream &os) const 
{

   os << "\nMBEulerMM::printClassData..." << endl;
   os << "MBEulerMM: this = " << (MBEulerMM*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::MultiblockDeformingGridGeometry<NDIM>*)d_grid_geometry << endl;

   os << "Parameters for physical problem ..." << endl;
   os << "   d_gamma = " << d_gamma << endl;
   os << "   d_alpha = " << d_alpha << endl;
   os << "   d_beta  = " << d_beta  << endl;
   os << "   d_xlo   = " << d_xlo[0] << ", " << d_xlo[1] << endl;
   os << "   d_xhi   = " << d_xhi[0] << ", " << d_xhi[1] << endl;
   os << "   d_ncells= " << d_ncells << endl;

}

/*
*************************************************************************
* 从输入数据库读取数据. 如果从输入数据库和重启动数据库中读取相同数据成  *
* 员的值.那么用输入数据库中的值覆盖重启动数据库的值.                    *
*************************************************************************
*/

void MBEulerMM::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if(!is_from_restart) {

      if (db->keyExists("gamma") ) {
         d_gamma = db->getDouble("gamma");
      }
      if (db->keyExists("alpha") ) {
         d_alpha = db->getDouble("alpha");
      }
      if (db->keyExists("beta") ) {
         d_beta = db->getDouble("beta");
      }
      if (db->keyExists("xlo") ) {
         db->getDoubleArray("xlo",d_xlo,NDIM);
      } else {
         TBOX_ERROR("no key `xlo' is found.");
      }
      if (db->keyExists("xhi") ) {
         db->getDoubleArray("xhi",d_xhi,NDIM);
      } else {
         TBOX_ERROR("no key `xhi' is found.");
      }
      if (db->keyExists("ncells") ) {
         db->getIntegerArray("ncells",d_ncells,NDIM);
      } else {
         TBOX_ERROR("no key `ncells' is found.");
      }
   }

}

/*
*************************************************************************
* 将MBEulerMM类的数据成员写入到重启动数据库.                                *
*************************************************************************
*/

void MBEulerMM::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putDouble("d_gamma", d_gamma);
   db->putDouble("d_alpha", d_alpha);
   db->putDouble("d_beta", d_beta);
   db->putDoubleArray("xlo",d_xlo,NDIM);
   db->putDoubleArray("xhi",d_xhi,NDIM);
   db->putIntegerArray("ncells",d_ncells,NDIM);

}

/*
*************************************************************************
*    从重启动数据库中读取数据.                                          *
*************************************************************************
*/
void MBEulerMM::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
   }

   d_gamma = db->getDouble("d_gamma");
   d_alpha = db->getDouble("d_alpha");
   d_beta  = db->getDouble("d_beta");
   db->getDoubleArray("xlo",d_xlo,NDIM);
   db->getDoubleArray("xhi",d_xhi,NDIM);
   db->getIntegerArray("ncells",d_ncells,NDIM);

}


#endif


