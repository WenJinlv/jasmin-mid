 //
// 文件名:  MBEulerLevelIntegratorMM.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  采用移动网格方法, 求解Euler方程组.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "MBEulerLevelIntegratorMM.h" 
#include "MBEulerMM.h" 

#include "MultiblockPatchLevel.h" 
#include "tbox/MPI.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"


using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
MBEulerLevelIntegratorMM::MBEulerLevelIntegratorMM(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        MBEulerMM* patch_strategy,
        tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
        const bool register_for_restart)
:
   d_object_name(object_name),
   d_registered_for_restart(register_for_restart),
   d_patch_strategy(patch_strategy),
   d_grid_geometry(grid_geometry)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!input_db.isNull());
   assert(!grid_geometry.isNull());
   assert(patch_strategy !=
           ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

   if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->
        registerRestartItem(d_object_name, this);
   }

   d_cfl = 0.5;

   bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if(from_restart) getFromRestart();
   getFromInput(input_db, from_restart);

}

/****************************************************************
* 析构函数.
*****************************************************************/
MBEulerLevelIntegratorMM::~MBEulerLevelIntegratorMM() {
    if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void MBEulerLevelIntegratorMM::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
	{

	  // 初值构件 : 初始化新创建的网格层.
	  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              true, manager);
	 
	  // 步长构件 : 计算时间步长.
	  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                              d_patch_strategy,
                                                              manager);

	  // 数值构件: 计算通量,守恒差分. 
	  d_grad_intc = new algs::NumericalIntegratorComponent<NDIM>("COMPUTE_GRAD",
                                                              d_patch_strategy,
                                                              manager);
	  d_pde_intc  = new algs::NumericalIntegratorComponent<NDIM>("SOLVE_EQUATION",
                                                              d_patch_strategy,
                                                              manager);

	  // 数值构件: 网格移动.
	  d_wc0_intc  = new algs::NumericalIntegratorComponent<NDIM>("COMPUTE_WCNEW",
                                                              d_patch_strategy,
                                                              manager);
	  d_wci_intc  = new algs::NumericalIntegratorComponent<NDIM>("WCNEW_ITERATION",
                                                              d_patch_strategy,
                                                              manager);
	  d_coords_intc = new algs::NumericalIntegratorComponent<NDIM>("MOVING_MESH",
                                                              d_patch_strategy,
                                                              manager);
	  d_movegrad_intc = new algs::NumericalIntegratorComponent<NDIM>("MOVING_GRAD",
                                                              d_patch_strategy,
                                                              manager);
	  d_movesolu_intc = new algs::NumericalIntegratorComponent<NDIM>("MOVING_SOLUTION",
                                                              d_patch_strategy,
                                                              manager);

	  // 复制构件: 接收数值解, 同步前复制数据片.
	  d_reset_intc   = new algs::CopyIntegratorComponent<NDIM>("RESET_SOLUTION",
                                                              d_patch_strategy,
                                                              manager);
	  d_cycle_copy_intc = new algs::CopyIntegratorComponent<NDIM>("COPY_PREMOVING",
                                                              d_patch_strategy,
                                                              manager);

	  // 内存构件: 当前值数据片, 新值数据片, 演算数据片.
	  d_new_intc     = new algs::MemoryIntegratorComponent<NDIM>("NEW_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);
	  d_scratch_intc = new algs::MemoryIntegratorComponent<NDIM>("SCRATCH_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);

          // 规约构件: 计算迭代误差.
          d_iter_error_intc = new algs::ReductionIntegratorComponent<NDIM>(
                                                             "ITER_ERROR",
                                                              MPI_MAX,
                                                              d_patch_strategy,
                                                              manager);

}

/****************************************************************
* 初始化网格层.
*****************************************************************/
void MBEulerLevelIntegratorMM::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     NULL_USE(can_be_refined);

     // 初始化网格层.
     d_init_intc->initializeLevelData(hierarchy,
                                      level_number,
                                      init_data_time,
                                      initial_time,
                                      old_level,
                                      false);

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double MBEulerLevelIntegratorMM::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
     double dt = d_dt_intc->getLevelDt(level,
                                       dt_time,
                                       initial_time,
                                       flag_last_dt,
                                       last_dt,
                                       false);

     return(dt*d_cfl);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int MBEulerLevelIntegratorMM::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(predict_dt>=min_dt && predict_dt<=max_dt);
#endif

    NULL_USE(hierarchy_step_number);

    // 显式时间离散格式, 取时间步长等于预测时间步长.
    actual_dt = predict_dt;

    // 重构后或时间步序列的第1步, 为新值数据片和演算数据片调度内存空间.
    if(first_step) {
       d_new_intc->allocatePatchData(level,current_time);
       d_scratch_intc->allocatePatchData(level,current_time);
    } else {
       d_new_intc->setTime(level,current_time);
       d_scratch_intc->setTime(level,current_time);
    }

    // 复制数据: current->new.
    d_cycle_copy_intc-> copyPatchData(level,
                                      current_time);

    // 移动网格 : 循环迭代.
    bool cycle = true;
    int ITS = 0;
    while(cycle) {

       // 移动网格: 第1步, 为控制函数wc赋初值.
       d_wc0_intc -> computing(level,
                               current_time,
                               actual_dt);

       // 移动网格: 第2步, 控制函数迭代.
       for(int i=0; i<4; i++) {

           d_wci_intc -> computing(level,
                                   current_time,
                                   actual_dt);
       }

       // 移动网格: 第3步, 移动网格结点.
       d_coords_intc -> computing(level,
                                  current_time,
                                  actual_dt); 

      // 移动网格: 第4步, 计算梯度.
      d_movegrad_intc->computing(level,
                                 current_time,
                                 actual_dt,
                                 false);

       // 移动网格: 第4步, 网格移动后, 修正数值解.
       d_movesolu_intc -> computing(level,
                                current_time,
                                actual_dt);

       // 移动网格: 第5步, 统计误差.
       double error = 0.0;
       d_iter_error_intc->reduction(&error,1,level,current_time,actual_dt,false);

       ITS++;
       if(error < 1.e-4 || ITS > 10) cycle = false; 

       #ifdef DEBUG_CHECK_ASSERTIONS
       tbox::pout << "After " << ITS << " iterations , "
                  << "error = " << error << endl;
       #endif

    }

    // 计算梯度.
    d_grad_intc->computing(level,
                           current_time,
                           actual_dt);

    // 解Euler方程.
    d_pde_intc->computing(level,
                          current_time,
                          actual_dt);
                  
    // 设置新值数据片的时刻.
    d_new_intc->setTime(level,current_time+actual_dt);

    return(1);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void MBEulerLevelIntegratorMM::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{
    // 将数值解从新值数据片复制到当前值数据片.
    d_reset_intc->copyPatchData(level,
                                new_time);

    // 释放新值数据片.
    if(deallocate_data) d_new_intc->deallocatePatchData(level);

}


/****************************************************************
* 打印. 
*****************************************************************/
void MBEulerLevelIntegratorMM::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
          << (MBEulerLevelIntegratorMM*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_cfl = " << d_cfl << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart <<endl;

    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::GridGeometry<NDIM>*)d_patch_strategy << endl;

}

/****************************************************************
* 从输入文件中读取参数值.
*****************************************************************/
void MBEulerLevelIntegratorMM::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    if (db->keyExists("cfl")) {
        d_cfl = db->getDouble("cfl");
    } else {
        if (!is_from_restart) {
            d_cfl = db->getDoubleWithDefault("cfl", d_cfl);
        }
    }

}

/****************************************************************
* 从重启动数据库中读取参数值.
*****************************************************************/
void MBEulerLevelIntegratorMM::getFromRestart()
{
    tbox::Pointer<tbox::Database> root_db =
        tbox::RestartManager::getManager()->getRootDatabase();

    tbox::Pointer<tbox::Database> db;
    if ( root_db->isDatabase(d_object_name) ) {
        db = root_db->getDatabase(d_object_name);
    } else {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file" << endl);
    }

    d_cfl = db->getDouble("d_cfl");

}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void MBEulerLevelIntegratorMM::putToDatabase(
         tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    db->putDouble("d_cfl", d_cfl);

}
