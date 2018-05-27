 //
// 文件名:  EulerLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  Euler方程组网格层时间积分算法的实现.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "EulerLevelIntegrator.h" 

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
EulerLevelIntegrator::EulerLevelIntegrator(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
        tbox::Pointer< hier::GridGeometry<NDIM> > grid_geometry,
        const bool register_for_restart,
        const bool use_time_refinement)
:
   d_object_name(object_name),
   d_use_time_refinement(use_time_refinement),
   d_registered_for_restart(register_for_restart),
   d_patch_strategy(patch_strategy),
   d_grid_geometry(grid_geometry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geometry.isNull());
   assert(patch_strategy != ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

   if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->
        registerRestartItem(d_object_name, this);
   }

   d_cfl = 0.9;

   bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if(from_restart) getFromRestart();
   getFromInput(input_db, from_restart);

}

/****************************************************************
* 析构函数.
*****************************************************************/
EulerLevelIntegrator::~EulerLevelIntegrator() {
    if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void EulerLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

  // 初值构件 : 初始化新创建的网格层.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              false, manager);
 
  // 步长构件 : 计算时间步长.
  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                        d_patch_strategy,
                                                        manager);

  // 数值构件: 计算通量,守恒差分,标记. 
  d_flux_intc = new algs::NumericalIntegratorComponent<NDIM>("FLUX",
                                                        d_patch_strategy,
                                                        manager);
  d_diff_intc = new algs::NumericalIntegratorComponent<NDIM>("CONSERVATIVE_DIFFERENCE",
                                                        d_patch_strategy,
                                                        manager);
  d_tag_intc  = new algs::NumericalIntegratorComponent<NDIM>("TAGC",
                                                        d_patch_strategy,
                                                        manager);

  // 同步构件: 粗细网格层交接面流连续, 细网格层覆盖的粗网格层内部区域.
  d_sync_flux_intc    = new algs::SynchronizeIntegratorComponent<NDIM>("SYNC_FLUX","INTERFACE_SYNC",
                                                              d_patch_strategy,
                                                              manager);
  d_sync_overlay_intc = new algs::SynchronizeIntegratorComponent<NDIM>("SYNC_OVERLAY","OVERLAY_SYNC",
                                                              d_patch_strategy,
                                                              manager);

  // 复制构件: 接收数值解, 同步前复制数据片.
  d_reset_intc         = new algs::CopyIntegratorComponent<NDIM>("RESET_SOLUTION",
                                                              d_patch_strategy,
                                                              manager);
  d_sync_copy_intc     = new algs::CopyIntegratorComponent<NDIM>("SYNC_COPY",
                                                              d_patch_strategy,
                                                              manager);

  // 内存构件: 当前值数据片, 新值数据片, 演算数据片.
  d_new_intc     = new algs::MemoryIntegratorComponent<NDIM>("NEW_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);
  d_scratch_intc = new algs::MemoryIntegratorComponent<NDIM>("SCRATCH_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);

}

/****************************************************************
* 初始化网格层.
*****************************************************************/
void EulerLevelIntegrator::initializeLevelData(
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
                                      old_level);

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double EulerLevelIntegrator::getLevelDt(
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
                                       last_dt);

     return(dt*d_cfl);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int EulerLevelIntegrator::advanceLevel(
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

    // 调用标准积分函数, 完成1个时间步的积分.
    standardAdvanceLevel(level,
                         current_time,
                         actual_dt,
                         first_step);

    return(1);

}

/****************************************************************
* 细网格层校正相邻粗网格层, 积分一个时间步长.
*****************************************************************/
void EulerLevelIntegrator::synchronizeAdvanceLevel(
           const tbox::Pointer< hier::BasePatchLevel<NDIM> > coarser_level, 
           const tbox::Pointer< hier::BasePatchLevel<NDIM> > finer_level, 
           const double current_time,
           const double actual_dt)
{

    // 校正通量.
    d_sync_flux_intc->synchronizeCoarserLevel(finer_level,
                                              current_time+actual_dt);

    // 复制当前值数据片到新值数据片.
    d_sync_copy_intc->copyPatchData(coarser_level,
                                    current_time+actual_dt);
    
    // 守恒差分.
    d_diff_intc->computing(coarser_level,
                           current_time,
                           actual_dt,
                           false);

}

/****************************************************************
* 时间步序列中, 积分一个时间步长.
*****************************************************************/
void EulerLevelIntegrator::standardAdvanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double actual_dt,
                const bool   first_step)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif

    // 重构后或时间步序列的第1步, 为新值数据片调度内存空间.
    if(first_step) {
       d_new_intc->allocatePatchData(level,
                                     current_time);
    } else {
       d_new_intc->setTime(level,current_time);
    }
 
    // 计算通量. 
    d_flux_intc->computing(level,
                           current_time,
                           actual_dt);

    // 存储通量, 以便细网格层校正相邻的粗网格层, 保证流连续连接.
    if(level->getLevelNumber()>0) {
       bool initialize_fluxsum = first_step && d_use_time_refinement ||
                                 !d_use_time_refinement;
       d_sync_flux_intc->storeInterfaceFluids(level,
                                              initialize_fluxsum);
    }

    // 守恒差分.
    d_diff_intc->computing(level,
                           current_time,
                           actual_dt,
                           false);

    // 设置新值数据片的时刻.
    d_new_intc->setTime(level,current_time+actual_dt);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void EulerLevelIntegrator::acceptTimeDependentSolution(
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
* 同步新的网格片层次结构.
*****************************************************************/
void  EulerLevelIntegrator::synchronizeNewHierarchy(
                   const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                   const int coarsest_level, 
                   const int finest_level, 
                   const double sync_time, 
                   const bool initial_time)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   NULL_USE(sync_time);
   NULL_USE(initial_time);
};


/****************************************************************
* 采用细化时间积分模式.
*****************************************************************/
bool EulerLevelIntegrator::usingRefinedTimestepping()
{
   return(d_use_time_refinement);
}

/****************************************************************
* 同步粗网格层.
*****************************************************************/
void EulerLevelIntegrator::synchronizeCoarserLevel(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int finer_level_number,
      const double sync_time, 
      const double coarser_old_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(finer_level_number>0);
#endif

   // 校正相邻粗网格层的通量. 
   tbox::Pointer< hier::PatchLevel<NDIM> > finer_level =
                 hierarchy -> getPatchLevel(finer_level_number);
   tbox::Pointer< hier::PatchLevel<NDIM> > coarser_level =
                 hierarchy -> getPatchLevel(finer_level_number-1);

   const double actual_dt = sync_time-coarser_old_time;
   synchronizeAdvanceLevel(coarser_level,
                           finer_level,
                           coarser_old_time,
                           actual_dt);

   // 校正细网格层覆盖的计算区域.
   d_sync_overlay_intc->synchronizeCoarserLevel(finer_level, sync_time);

}

/****************************************************************
* 打标记.
*****************************************************************/
void EulerLevelIntegrator::tagCellsForRefinement(
               const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
               const int level_number,
               const double error_data_time,
               const int tag_index,
               const bool initial_time)
{

   d_scratch_intc->allocatePatchData(
                                  hierarchy->getPatchLevel(level_number),
                                  error_data_time);

   // 设置存储标记值的数据片索引号.
   d_patch_strategy->setTagIndex(tag_index);
   d_tag_intc->computing(hierarchy->getPatchLevel(level_number),
                         error_data_time,
                         0.0, true, initial_time);
   d_patch_strategy->clearTagIndex();

   d_scratch_intc->deallocatePatchData(
                         hierarchy->getPatchLevel(level_number));

}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void EulerLevelIntegrator::putToDatabase(
         tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    db->putDouble("d_cfl", d_cfl);

}

/****************************************************************
* 打印. 
*****************************************************************/
void EulerLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
          << (EulerLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_cfl = " << d_cfl << "\n"
       << "d_use_time_refinement = " << d_use_time_refinement << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart <<endl;

    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::GridGeometry<NDIM>*)d_patch_strategy << endl;

    os << "Integrator Components listed : \n" 
       << "  Initialize Integrator Component = " 
             << d_init_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_init_intc << "\n"
       << "  Step size Integrator Component = " 
             << d_dt_intc->getName() << " " 
             << (algs::IntegratorComponent<NDIM> *)d_dt_intc << "\n"
       << "  Numerical Integrator Component for Flux = " 
             << d_flux_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_flux_intc << "\n"  
       << "  Numerical Integrator Component for Cons.Diff. = " 
             << d_diff_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_diff_intc << "\n"  
       << "  Tag Integrator Component = " 
             << d_tag_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_tag_intc << "\n"  
       << "  Synchronize Integrator Component for Interface = " 
             << d_sync_flux_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_sync_flux_intc << "\n"
       << "  Synchronize Integrator Component for Overlay = " 
             << d_sync_overlay_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_sync_overlay_intc << "\n"  
       << "  Copy Integrator Component for accept solution = " 
             << d_reset_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_reset_intc << "\n"
       << "  Memory Integrator Component for new solution = " 
             << d_new_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_new_intc << "\n"
       << "  Memory Integrator Component for scratch = " 
             << d_scratch_intc->getName() << " "
             << (algs::IntegratorComponent<NDIM> *)d_scratch_intc << "\n"  
      << endl;
}

/****************************************************************
* 从输入文件中读取参数值.
*****************************************************************/
void EulerLevelIntegrator::getFromInput(
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
void EulerLevelIntegrator::getFromRestart()
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
