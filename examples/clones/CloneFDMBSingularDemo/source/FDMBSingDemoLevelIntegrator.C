 //
// 文件名:  FDMBSingDemoLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  实现网格层时间积分算法.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "FDMBSingDemoLevelIntegrator.h" 
#include "FDMBSingDemo.h" 

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
FDMBSingDemoLevelIntegrator::FDMBSingDemoLevelIntegrator(
        const string& object_name,
        algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
        tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry)
:
   d_object_name(object_name),
   d_patch_strategy(patch_strategy),
   d_grid_geometry(grid_geometry)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(patch_strategy !=
           ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

}

/****************************************************************
* 析构函数.
*****************************************************************/
FDMBSingDemoLevelIntegrator::~FDMBSingDemoLevelIntegrator() {
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void FDMBSingDemoLevelIntegrator::initializeLevelIntegrator(
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
 
  // 数值构件 : 完成累加操作.
  d_sum_intc   = new algs::NumericalIntegratorComponent<NDIM>("SUMMING",
                                                              d_patch_strategy,
                                                              manager);
 
  // 复制构件 : 接收数值解.
  d_reset_intc   = new algs::CopyIntegratorComponent<NDIM>("ACCEPT_SOLUTION",
                                                              d_patch_strategy,
                                                              manager);

  // 内存构件: 新值数据片.
  d_new_intc  = new algs::MemoryIntegratorComponent<NDIM>("NEW_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);

  // 内存构件: 演算数据片.
  d_alloc_scratch_data = new algs::MemoryIntegratorComponent<NDIM>(
                                           "ALLOC_SCRATCH_DATA",
                                           d_patch_strategy,
                                           manager);

  // 克隆广播构件:将数值解从隆主网格层广播到克隆网格层.
  d_bcast_solution = new algs::CloneBcastComponent<NDIM>(
                                                       "BCAST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  // 克隆汇集构件:将数值解从克隆网格层汇集到属主网格层.
  d_reduction_solution = new algs::CloneReductionComponent<NDIM>(
                                                       "REDUCTION_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  // 数值构件: 前处理属主网格层.
  d_prepprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "PREP_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

  // 数值构件: 后处理属主网格层.
  d_postprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "POST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);


}

/****************************************************************
* 初始化网格层.
*****************************************************************/
void FDMBSingDemoLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     NULL_USE(can_be_refined);

     // 初始化属主网格层.
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
double FDMBSingDemoLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{

//???
     double dt;
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
int FDMBSingDemoLevelIntegrator::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double predict_dt,
                const bool   first_step)
{
    hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
    {
        // 重构后或时间步序列的第1步, 为新值数据片和演算数据片调度内存空间.
        if(first_step) {
          d_new_intc->allocatePatchData(level,current_time);
        } else {
           d_new_intc->setTime(level,current_time);
        }
    }
    hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

    // 前处理属主网格层.
    d_prepprocess_computation->computing(level,current_time,0);

    hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
    {
        // 广播数值解.
        d_bcast_solution->bcast(level);

        // 数值构件.
        d_sum_intc->computing(level,
                              current_time,
                              predict_dt);
    
        // 设置新值数据片的时刻.
        d_new_intc->setTime(level,current_time+predict_dt);

        // 汇集数值解.
        d_reduction_solution->reduce(level, xfer::RefineUtilities::SUM);

     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

    return(1);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void FDMBSingDemoLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{
    // 将数值解从新值数据片复制到当前值数据片.
    d_reset_intc->copyPatchData(level,
                                new_time);

    // 后处理数值解.
    d_postprocess_computation->computing(level,new_time,0);

    // 释放新值数据片.
    hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
    if(deallocate_data) d_new_intc->deallocatePatchData(level);
    hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();
}

/*************************************************************************
 * 为邦员之间的数据传递调度内存空间.
 ************************************************************************/
void FDMBSingDemoLevelIntegrator::allocateFederalXferData(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!level.isNull());
#endif

     // 开辟演算数据片的内存.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     d_alloc_scratch_data->allocatePatchData(level, time);
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

}

/*************************************************************************
 * 为邦员之间的数据传递调度内存空间.
 ************************************************************************/
void FDMBSingDemoLevelIntegrator::deallocateFederalXferData(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!level.isNull());
#endif

     // 释放演算数据片的内存.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     d_alloc_scratch_data->deallocatePatchData(level);
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

}


/****************************************************************
* 打印. 
*****************************************************************/
void FDMBSingDemoLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
          << (FDMBSingDemoLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;

    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::MultiblockGridGeometry<NDIM>*)d_grid_geometry<< endl;

}

