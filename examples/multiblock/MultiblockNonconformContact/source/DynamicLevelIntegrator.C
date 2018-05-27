//
// 文件名:  DynamicLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  类DynamicLevelIntegrator的实现.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include <assert.h>

#include "Dynamic.h"
#include "DynamicLevelIntegrator.h"

#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Tracer.h"
#include "tbox/Statistician.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

// 声明Fortran数值子程序的头文件.
#include "DynamicFort.h"

using namespace JASMIN;

/****************************************************************
 * 构造函数.                                                    *
 ****************************************************************/
DynamicLevelIntegrator::DynamicLevelIntegrator(
  const string& object_name,
  tbox::Pointer<tbox::Database> input_db,
  Dynamic* patch_strategy,
  tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
  tbox::Pointer< hier::MultiblockPatchHierarchy<NDIM> > patch_hierarchy,
  const bool register_for_restart,
  const bool use_time_refinement,
  const bool use_implicit_discretization)
    :
    d_object_name(object_name),
    d_registered_for_restart(register_for_restart),
    d_patch_strategy(patch_strategy),
    d_grid_geometry(grid_geometry),
    d_patch_hierarchy(patch_hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::DynamicLevelIntegrator");
  assert(!object_name.empty());
  assert(!input_db.isNull());
  assert(!grid_geometry.isNull());
  assert(!patch_hierarchy.isNull());
  assert(patch_strategy != ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

  NULL_USE(use_implicit_discretization);

  if (d_registered_for_restart)
    {
      tbox::RestartManager::getManager()->
      registerRestartItem(d_object_name, this);
    }

  bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
  if (from_restart)
    {
      getFromRestart();
    }
  getFromInput(input_db, from_restart);

  t_nonconform_contact_computing = tbox::TimerManager::getManager()->
        getTimer("apps::DynamicLevelIntegrator::nonconform_contact_computing");
}

/****************************************************************
 * 析构函数.                                                    *
 ****************************************************************/
DynamicLevelIntegrator::~DynamicLevelIntegrator()
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::~DynamicLevelIntegrator");
#endif

  if (d_registered_for_restart)
    {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
}

/****************************************************************
 * 创建所有积分构件.                                            *
 ****************************************************************/
void DynamicLevelIntegrator::initializeLevelIntegrator(
  tbox::Pointer<algs::IntegratorComponentManager<NDIM> > component_manager)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::initializeLevelIntegrator");
#endif

  // 创建初始化构件, 用于初始化新创建的网格层.
  d_init_component = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                     d_patch_strategy,
                     true,
                     component_manager);

  // 创建步长构件, 用于计算时间步长.
  d_dt_component = new algs::DtIntegratorComponent<NDIM>("TIME_STEP",
                   d_patch_strategy,
                   component_manager);

  // 创建数值构件, 用于更新网格结点位置坐标.
  d_update_node_coordinate_component =
    new algs::NumericalIntegratorComponent<NDIM>("UPDATE_NODE_COORDINATE",
        d_patch_strategy,
        component_manager);

  // 创建非协调拼接构件, 用于验证非协调拼接构件正确性.
  d_nonconform_component =
    new algs::NonconformContactIntegratorComponent<NDIM>("NONCONFORM_CONTACT",
        d_patch_strategy,
        component_manager);

  // 创建复制构件, 用于将数据片从NONCONFORM上下文复制到CURRENT上下文
  d_accept_component = new algs::CopyIntegratorComponent<NDIM>("ACCEPT_SOLUTION",
                       d_patch_strategy,
                       component_manager);

   // 创建内存构件, 用于为NONCONFORM上下文中的数据片管理内存空间.
  d_memory_component = new algs::MemoryIntegratorComponent<NDIM>("NONCONFORM_PATCH_DATA",
                    d_patch_strategy,
                    component_manager);
}

/****************************************************************
 * 基于初始化构件, 初始化网格层数据.                            *
 ****************************************************************/
void DynamicLevelIntegrator::initializeLevelData(
  const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
  const int    level_number,
  const double init_data_time,
  const bool   can_be_refined,
  const bool   initial_time,
  const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
  const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::initializeLevelData");
#endif

  NULL_USE(can_be_refined);

  // 初始化数据片.
  d_init_component->initializeLevelData(hierarchy,
                                        level_number,
                                        init_data_time,
                                        initial_time,
                                        old_level,
                                        false);
}

/****************************************************************
 * 基于步长构件, 计算时间步长.                                  *
 ****************************************************************/
double DynamicLevelIntegrator::getLevelDt(
  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
  const double dt_time,
  const bool   initial_time,
  const int    flag_last_dt,
  const double last_dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::getLevelDt");
#endif

  // 时间步长在CURRENT上下文中计算.
  double dt = d_dt_component->getLevelDt(level,
                                         dt_time,
                                         initial_time,
                                         flag_last_dt,
                                         last_dt);
  //dt = dt * d_cfl;
  //tbox::pout << "d_cfl: " << d_cfl << "  dt: " << dt << endl;
  //double dt_grow = 1.1 * last_dt;
  //dt = tbox::Utilities::dmin(dt, dt_grow);
  //tbox::pout << "dt: " << dt << endl;

  return(dt);

}

/****************************************************************
 * 调用正常的网格层积分函数, 积分一个时间步长. 具体地, 根据速度 *
 * 更新网格结点位置坐标.                                        *
 ****************************************************************/
int DynamicLevelIntegrator::advanceLevel(
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
  tbox::Tracer tracer("DynamicLevelIntegrator::advanceLevel");
  assert((predict_dt >= min_dt) && (predict_dt <= max_dt));
#endif

  NULL_USE(hierarchy_step_number);
  // 显式时间离散格式, 取时间步长等于预测时间步长.
  actual_dt = predict_dt;

  if (first_step)
    {
      // 重构后或时间步序列的第一步, 为NONCONFORM上下文中的数据片调度内存空间.
      d_memory_component->allocatePatchData(level, current_time);
    }
  else
    {
      // 为NONCONFORM上下文中的数据片设置时戳.
      d_memory_component->setTime(level, current_time);
    }

#if 1
  // 调用数值构件, 更新网格结点坐标.
  d_update_node_coordinate_component->computing(level,
                                                current_time,
                                                actual_dt,
                                                false);
#endif
  // 非协调接触计算.
  t_nonconform_contact_computing->start();
  d_nonconform_component->computing(level,
                                   current_time,
                                   actual_dt,
                                   false,        // 不需要填充物理边界.
                                   first_step);

  t_nonconform_contact_computing->stop();

  // 设置新值数据片的时刻.
  d_memory_component->setTime(level,current_time+actual_dt);
  return(1);
}


/****************************************************************
 * 接收数值解.                                                  *
 ****************************************************************/
void DynamicLevelIntegrator::acceptTimeDependentSolution(
  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
  const double new_time,
  const bool deallocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::acceptTimeDependentSolution");
  assert(!level.isNull());
#endif

  // 将数值解从新值数据片(NEW上下文)复制到当前值数据片(CURRENT上下文).
  d_accept_component->copyPatchData(level,
                                    new_time);

  // 释放NONCONFORM上下文中数据片的内存空间.
  if (deallocate_data)
    {
      d_memory_component->deallocatePatchData(level);
    }
}

/****************************************************************
 * 输出数据成员到指定输出流.                                    *
 ****************************************************************/
void DynamicLevelIntegrator::printClassData(ostream& os) const
  {
    os << "\nDynamicLevelIntegrator::printClassData..." << endl;
    os << "DynamicTimeLevelIntegrator<NDIM>: this = "
    << (DynamicLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_registered_for_restart = " << d_registered_for_restart <<endl;
    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::GridGeometry<NDIM>*)d_patch_strategy << endl;
  }

/****************************************************************
 * 从输入文件中读取参数值.                                      *
 ****************************************************************/
void DynamicLevelIntegrator::getFromInput(
  tbox::Pointer<tbox::Database> db,
  bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("DynamicLevelIntegrator::getFromInput");
  assert(!db.isNull());
#endif
  NULL_USE (db);
  NULL_USE (is_from_restart);

}

/****************************************************************
 * 从重启动数据库中读取参数值.                                  *
 ****************************************************************/
void DynamicLevelIntegrator::getFromRestart()
{
  tbox::Pointer<tbox::Database> root_db =
    tbox::RestartManager::getManager()->getRootDatabase();

  tbox::Pointer<tbox::Database> db;
  if ( root_db->isDatabase(d_object_name) )
    {
      db = root_db->getDatabase(d_object_name);
    }
  else
    {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file" << endl);
    }

}

/****************************************************************
 * 输出数据成员到重启动数据库.                                  *
 ****************************************************************/
void DynamicLevelIntegrator::putToDatabase(
  tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif
}
