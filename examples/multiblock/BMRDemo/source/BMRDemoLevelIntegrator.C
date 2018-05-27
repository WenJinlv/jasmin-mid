 //
// 文件名:  BMRDemoLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  实现网格层时间积分算法.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "BMRDemoLevelIntegrator.h" 
#include "BMRDemo.h" 

#include "MultiblockPatchLevel.h" 
#include "tbox/MPI.h" 

#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"


using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
BMRDemoLevelIntegrator::BMRDemoLevelIntegrator(
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

   d_hierarchy_step_number = 0;
   d_current_time = 0.0;

}

/****************************************************************
* 析构函数.
*****************************************************************/
BMRDemoLevelIntegrator::~BMRDemoLevelIntegrator() {
    if (d_registered_for_restart) {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void BMRDemoLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

  // 初值构件 : 初始化新创建的网格层.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              true, // use_physical_domain_boxes
                                                              //false, // use_physical_domain_boxes
                                                              manager);
 
  // 步长构件 : 计算时间步长.
  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                              d_patch_strategy,
                                                              manager);
 
  // 数值构件 : 完成累加操作.
  d_comput_temperature   = new algs::NumericalIntegratorComponent<NDIM>("COMPUT_TEMP",
                                                              d_patch_strategy,
                                                              manager);
 
  // 复制构件 : 接收数值解.
  d_reset_intc   = new algs::CopyIntegratorComponent<NDIM>("ACCEPT_SOLUTION",
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
void BMRDemoLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     NULL_USE(can_be_refined);

    tbox::pout << "initializeLevelData in BMRDemoLevelIntegrator::initializeLevelData: begining..." << endl;
    tbox::pout << "  level_number in BMRDemoLevelIntegrator::initializeLevelData= " << level_number << endl;

    // 初始化网格层.
    d_init_intc->initializeLevelData(hierarchy,
                                     level_number,
                                     init_data_time,
                                     initial_time,
                                     old_level,
                                     false);

    tbox::pout << "initializeLevelData in BMRDemoLevelIntegrator::initializeLevelData: finished." << endl;

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double BMRDemoLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{

     tbox::pout << "BMRDemoLevelIntegrator::getLevelDt: begin." << endl;
     double dt = d_dt_intc->getLevelDt(level,
                                       dt_time,
                                       initial_time,
                                       flag_last_dt,
                                       last_dt,
                                       false);

     tbox::pout << "  dt = " << dt << endl;
     return(dt);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int BMRDemoLevelIntegrator::advanceLevel(
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

    tbox::pout << "BMRDemoLevelIntegrator::advanceLevel =" << first_step << endl;

    NULL_USE(hierarchy_step_number);
    d_hierarchy_step_number = hierarchy_step_number;
    d_current_time = current_time;

    // 显式时间离散格式, 取时间步长等于预测时间步长.
    actual_dt = predict_dt;

    tbox::pout << "  first_step =" << first_step << endl;
    // 重构后或时间步序列的第1步, 为新值数据片和演算数据片调度内存空间.
    if(first_step) {
       d_new_intc->allocatePatchData(level,current_time);
    } else {
       d_new_intc->setTime(level,current_time);
    }

    d_scratch_intc->allocatePatchData(level,current_time);
    tbox::pout << "  d_comput_temperature..." << endl;
    // 数值构件.
    d_comput_temperature->computing(level,
                          current_time,
                          actual_dt);

    // 设置新值数据片的时刻.
    d_new_intc->setTime(level,current_time+actual_dt);

    d_scratch_intc->deallocatePatchData(level);

    return(1);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void BMRDemoLevelIntegrator::acceptTimeDependentSolution(
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
void BMRDemoLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
          << (BMRDemoLevelIntegrator*)this << endl;
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
void BMRDemoLevelIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{
}

/****************************************************************
* 从重启动数据库中读取参数值.
*****************************************************************/
void BMRDemoLevelIntegrator::getFromRestart()
{
}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void BMRDemoLevelIntegrator::putToDatabase(
         tbox::Pointer<tbox::Database> db)
{
}


/*************************************************************************
 *
 * 标记待细化的网格块.
 *
 ************************************************************************/
void BMRDemoLevelIntegrator::tagBlocksForBMR(
               const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
               tbox::Array<int> tag_blocks_indx )
{

   tbox::pout << "BMRDemoLevelIntegrator::tagBlocksForBMR...." << endl;

   int nblocks = tag_blocks_indx.size();
   tbox::pout << "  d_current_time= " << d_current_time << endl;
   tbox::pout << "  d_hierarchy_step_number= " << d_hierarchy_step_number << endl;
   tbox::pout << "  mod(d_hierarchy_step_number, nblocks)= " << d_hierarchy_step_number % nblocks << endl;

   /* 为每块网格打标记是否需要细化或粗化. */
   for (int nb=0; nb<nblocks; nb++) {

#if 0

     if (d_hierarchy_step_number <= 2*nblocks-1) {
       if (nb==d_hierarchy_step_number % nblocks) {
          tag_blocks_indx[nb] = 1;
       } else if (d_hierarchy_step_number <= nblocks) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }
     } else {
       if (nb==d_hierarchy_step_number % nblocks) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }
     }

#endif

#if 1
     if (d_current_time >= 150.0 && d_current_time < 400.0) {

       if (nb == 1) {
          tag_blocks_indx[nb] = 1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 400.0 && d_current_time < 700.0) {

       if (nb == 2 || nb == 3) {
          tag_blocks_indx[nb] = 1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 700.0 && d_current_time < 900.0) {

       if (nb == 4 || nb == 5) {
          tag_blocks_indx[nb] = 1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 900.0 && d_current_time < 1000.0) {

       if (nb == 6) {
          tag_blocks_indx[nb] = 1;
       } else if (nb == 0) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 1000.0 && d_current_time < 1500.0) {

       if (nb == 7 || nb == 8) {
          tag_blocks_indx[nb] = 1;
       } else if (nb == 6) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 1500.0 && d_current_time < 1600.0) {

       if (nb == 2 || nb == 3) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 1600.0 && d_current_time < 1700.0) {

       if (nb == 2 || nb == 3) {
          tag_blocks_indx[nb] = 1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else if (d_current_time >= 1700.0) {

       if (nb == 1) {
          tag_blocks_indx[nb] = -1;
       } else {
          tag_blocks_indx[nb] = 0;
       }

     } else {

       tag_blocks_indx[nb] = 0;

     }
#endif

   }

}

/*************************************************************************
 *
 * 细网格层同步相邻的粗网格层.
 *
 * 该函数用细网格层数据矫正粗网格层的守恒量 u. 
 ************************************************************************/
void BMRDemoLevelIntegrator::synchronizeCoarserLevel(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
      const int finer_level_number, 
      const double sync_time,
      const double coarser_old_time)
{ 
}

/*************************************************************************
 *
 * 标记待细化的网格单元。
 *
 ************************************************************************/
void BMRDemoLevelIntegrator::tagCellsForRefinement(
               const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
               const int level_number,
               const double error_data_time,
               const int tag_index,
               const bool initial_time)
{

}

/*************************************************************************
 *
 * 返回时间积分的步进模式.
 *
 ************************************************************************/
bool BMRDemoLevelIntegrator::usingRefinedTimestepping()
{
   return(false);
}

