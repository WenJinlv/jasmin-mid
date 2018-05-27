 //
// 文件名:  MBSingDemoLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  实现网格层时间积分算法.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "MBSingDemoLevelIntegrator.h" 
#include "MBSingDemo.h" 

#include "MultiblockPatchLevel.h" 
#include "tbox/MPI.h" 

#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "CommunicatorDatabase.h"

using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
MBSingDemoLevelIntegrator::MBSingDemoLevelIntegrator(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
        tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
        const bool register_for_restart)
:
   d_object_name(object_name),
   d_registered_for_restart(register_for_restart),
   d_patch_strategy(patch_strategy),
   d_grid_geometry(grid_geometry)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(patch_strategy !=
           ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

   if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->
        registerRestartItem(d_object_name, this);
   }

   bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if(from_restart) getFromRestart();
   getFromInput(input_db, from_restart);

}

/****************************************************************
* 析构函数.
*****************************************************************/
MBSingDemoLevelIntegrator::~MBSingDemoLevelIntegrator() {
    if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void MBSingDemoLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

  // 初值构件 : 初始化属主网格层.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              true, manager);

  //数值构件: 后处理属主网格层.
  d_postprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "POST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  // 步长构件 : 克隆网格层 计算时间步长.
  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                 d_patch_strategy,
                                                 manager);

  // 数值构件 : 克隆网格层完成累加操作.
  d_sum_intc   = new algs::NumericalIntegratorComponent<NDIM>("SUMMING",
                                                         d_patch_strategy,
                                                         manager);
 
  // 复制构件 : 克隆网格层接收数值解.
  d_reset_intc   = new algs::CopyIntegratorComponent<NDIM>("ACCEPT_SOLUTION",
                                                      d_patch_strategy,
                                                      manager);

  // 内存构件: 当前值数据片, 新值数据片, 演算数据片.
  d_current_intc = new algs::MemoryIntegratorComponent<NDIM>("CURRENT_ALLOC_PATCH_DATA",
                                                        d_patch_strategy,
                                                        manager);
  d_new_intc     = new algs::MemoryIntegratorComponent<NDIM>("NEW_ALLOC_PATCH_DATA",
                                                        d_patch_strategy,
                                                        manager);
  d_scratch_intc = new algs::MemoryIntegratorComponent<NDIM>("SCRATCH_ALLOC_PATCH_DATA",
                                                        d_patch_strategy,
                                                        manager);

  //克隆广播构件:将数值解从隆主网格层广播到克隆网格层.
  d_bcast_solution = new algs::CloneBcastComponent<NDIM>(
                                                       "BCAST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  //克隆汇集构件:将数值解从克隆网格层汇集到属主网格层.
  d_reduction_solution = new algs::CloneReductionComponent<NDIM>(
                                                       "REDUCTION_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  //数值构件: 前处理克隆网格层的数据.
  d_prepprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "PREP_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  // 外表面结点累加和构件
  d_onsum_intc   = 
    new algs::OuterdataOperationIntegratorComponent<NDIM>("OUTERNODE_SUMMING",
                                                     d_patch_strategy,
                                                     manager,
                                                     "SUM");

}

/****************************************************************
* 初始化网格层.
*****************************************************************/
void MBSingDemoLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     NULL_USE(can_be_refined);

     // 克隆计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {
         // 初始化属主网格层.
         d_init_intc->initializeLevelData(hierarchy,
                                      level_number,
                                      init_data_time,
                                      initial_time,
                                      old_level,
                                      false);

         // 克隆网格层调度内存空间, 用于存储属主网格层广播来的数值解. 
         d_current_intc->allocatePatchData(
               hierarchy->getPatchLevel(level_number),init_data_time);

         // 广播数值解.
         d_bcast_solution->init_bcast(hierarchy->getPatchLevel(level_number));

     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double MBSingDemoLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
     double dt;

     // 克隆计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {
        dt = d_dt_intc->getLevelDt(level,
                                       dt_time,
                                       initial_time,
                                       flag_last_dt,
                                       last_dt,
                                       false);
     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

     return(dt);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int MBSingDemoLevelIntegrator::advanceLevel(
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
    TBOX_ASSERT(predict_dt>=min_dt && predict_dt<=max_dt);
#endif

    NULL_USE(hierarchy_step_number);

    // 克隆计算.
    hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
    {
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

       // 广播数值解.
       d_bcast_solution->bcast(level);

       // 前处理数值解.
       d_prepprocess_computation->computing(level, current_time, predict_dt);

       // 数值构件.
       d_sum_intc->computing(level,
                             current_time,
                             actual_dt);

       // 外表面结点累加和.
       d_onsum_intc->operate(level);

       // 设置新值数据片的时刻.
       d_new_intc->setTime(level,current_time+actual_dt);

       // 汇集数值解.
       d_reduction_solution->reduce(level, xfer::RefineUtilities::SUM);

    }
    hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

    return(1);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void MBSingDemoLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{
    // 将数值解从新值数据片复制到当前值数据片.
    d_reset_intc->copyPatchData(level,
                                new_time);

    // 后处理数值解.
    d_postprocess_computation->computing(level,new_time,0);

    // 克隆计算: 释放新值数据片.
    hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
    if(deallocate_data) d_new_intc->deallocatePatchData(level);
    hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();
}

/****************************************************************
* 打印. 
*****************************************************************/
void MBSingDemoLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
          << (MBSingDemoLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_registered_for_restart = " << d_registered_for_restart <<endl;

    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::MultiblockGridGeometry<NDIM>*)d_grid_geometry<< endl;

}

/****************************************************************
* 从输入文件中读取参数值.
*****************************************************************/
void MBSingDemoLevelIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{
}

/****************************************************************
* 从重启动数据库中读取参数值.
*****************************************************************/
void MBSingDemoLevelIntegrator::getFromRestart()
{
}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void MBSingDemoLevelIntegrator::putToDatabase(
         tbox::Pointer<tbox::Database> db)
{
}
