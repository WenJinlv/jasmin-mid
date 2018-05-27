//
// 文件  :	DynamicsLevelIntegrator.C
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 1.20 $
// 修改  :	$Date: 2006/08/09 06:22:58 $
// 描述  :	网格层流体力学时间积分算法类的具体实现.
//

#include "DynamicsLevelIntegrator.h"

#include "VariableDatabase.h"
#include "tbox/IEEE.h"
#include "tbox/MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

/*
*************************************************************************
*                                                                       *
* 构造函数,初始化对象的数据成员.
*                                                                       *
*************************************************************************
*/

DynamicsLevelIntegrator::DynamicsLevelIntegrator(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   Dynamics* patch_strategy,
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
   bool register_for_restart) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(patch_strategy != ((Dynamics*)NULL));
#endif

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }
   
   d_patch_strategy = patch_strategy;
   
   d_grid_geometry = grid_geometry;

   d_cfl = 1.0;

   d_initial_dt = 1.0e-7;

   /*
    * 从输入数据库或重启动数据库中读取数据,初始化数据成员.
    */

   bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, from_restart);
}

/*
*************************************************************************
*                                                                       *
* 析构函数. 从重启动数据库中注销对象.
*                                                                       *
*************************************************************************
*/
DynamicsLevelIntegrator::~DynamicsLevelIntegrator()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void DynamicsLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

  // 初值构件 : 初始化新创建的网格层数据.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              true, manager);

  // 步长构件 : 计算时间步长.
  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                              d_patch_strategy,
                                                              manager);

  // 数值构件 : 计算压力和黏性
  d_preprocess_intc = new algs::NumericalIntegratorComponent<NDIM>("PRE",
                                                              d_patch_strategy,
                                                              manager);

  // 数值构件 : 计算网格边上的压力
  d_side_pressure_intc = new algs::NumericalIntegratorComponent<NDIM>("SIDE_PRESSURE",
                                                              d_patch_strategy,
                                                              manager);

  // 数值构件 : 计算加速度和速度
  d_velocity_intc = new algs::NumericalIntegratorComponent<NDIM>("VELOCITY",
                                                              d_patch_strategy,
                                                              manager);

  // 数值构件 : 计算位置、密度和能量
  d_grid_dens_energy_intc = new algs::NumericalIntegratorComponent<NDIM>("ENERGY",
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

  // 外表面结点量舍入误差构件.
  d_outernode_roundoff = 
     new algs::OuterdataOperationIntegratorComponent<NDIM>("ROUNDOFF",
                                                           d_patch_strategy,
                                                           manager,
                                                           "ROUND_OFF");

}



/*
*************************************************************************
*                                                                       *
* 在网格层上,初始化所有数据片.
*
* 该函数按如下步骤执行:  
*     (1) 为指定网格层上的数据片调度内存.
*     (2) 利用初值构件, 初始化数据片数据.
*
*************************************************************************
*/

void DynamicsLevelIntegrator::initializeLevelData(
   const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
   const int    level_number,
   const double init_data_time,
  const bool   can_be_refined,
  const bool   initial_time,
  const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
  const bool allocate_data)
{
     // 初始化网格层.
     d_init_intc->initializeLevelData(hierarchy,
                                      level_number,
                                      init_data_time,
                                      initial_time,
                                      old_level,
                                      false);

}

/*
*************************************************************************
*                                                                       *
* 返回网格层允许的最大时间步长.
* 该函数利用步长构件实现.
*
*************************************************************************
*/

double DynamicsLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
              const double dt_time,
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
   double dt;
   if(initial_time){
     dt = d_initial_dt; 
   }else{
     dt  = d_dt_intc->getLevelDt(level,dt_time,initial_time,flag_last_dt,last_dt,false);
     dt *= d_cfl;
   }
   
   return(dt);
}


/*
*************************************************************************
*                                                                       *
* 在指定的网格层上,使用当前的时间步长,积分守恒方程,
* 将流体力学数值解从当前时刻积分到新的时刻.
* 
* 该函数按如下步骤执行:  
* (1) 调用数值构件, 预处理本次积分所需的辅助状态量(压力和黏性).
* (2) 调用数值构件, 网格边上的压力. 
* (3) 调用数值构件, 计算网格结点的加速度(IGA格式)和速度.
* (4) 调用舍入误差构件, 消除网格片外表面结点量的差异.
* (5) 调用数值构件, 计算网格结点坐标、网格单元密度和能量.
*                                                                       *
*************************************************************************
*/

int DynamicsLevelIntegrator::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt)
{
    // 重构后的第1步, 为新值数据片和演算数据片调度内存空间.
    if(first_step) {
       d_new_intc->allocatePatchData(level,current_time);
       d_scratch_intc->allocatePatchData(level,current_time);
    } else {
       d_new_intc->setTime(level,current_time);
       d_scratch_intc->setTime(level,current_time);
    }

    // 显式时间离散格式, 取时间步长等于预测时间步长.
    actual_dt = predict_dt;


    // 数值构件.
    d_preprocess_intc->computing(level,
                          current_time,
                          actual_dt);

    // 数值构件.计算边压力
    d_side_pressure_intc->computing(level,
                          current_time,
                          actual_dt);

    // 数值构件.
    d_velocity_intc->computing(level,
                          current_time,
                          actual_dt);

#if ERROR_ROUNDOFF
    // 在舍入误差的可控范围之内, 消除网格片外表面结点量的差异.
    d_outernode_roundoff->operate(level);
#endif

    d_grid_dens_energy_intc->computing(level,
                          current_time,
                          actual_dt);


    // 设置新值数据片的时刻.
    d_new_intc->setTime(level,current_time+actual_dt);

    
    return(1);
}

/*
*************************************************************************
*                                                                       *
* 接收数值解.
* 该函数利用复制构件实现.
*                                                                       *
*************************************************************************
*/

void DynamicsLevelIntegrator::acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                  const double new_time,
                  const bool deallocate_data)
{
    // 将数值解从新值数据片复制到当前值数据片.
    d_reset_intc->copyPatchData(level, new_time);

}



/*
*************************************************************************
*                                                                       *
* 从输入数据库中读取参数值.						*
*                                                                       *
*************************************************************************
*/

void DynamicsLevelIntegrator::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("cfl")) {
      d_cfl = db->getDouble("cfl");
   }

   if (db->keyExists("initial_dt")) {
      d_initial_dt = db->getDouble("initial_dt");
   }

}

/*
*************************************************************************
*                                                                       *
* 将数据成员写入重启动数据库.
*                                                                       *
*************************************************************************
*/

void DynamicsLevelIntegrator::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putDouble("d_cfl", d_cfl);
   db->putDouble("d_initial_dt", d_initial_dt);

}


/*
*************************************************************************
*                                                                       *
* 从重启动数据库中读取数据,初始化数据成员.
*                                                                       *
*************************************************************************
*/
void DynamicsLevelIntegrator::getFromRestart()
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
   d_initial_dt= db->getDouble("d_initial_dt");

}

/*
*************************************************************************
*                                                                       *
* 打印输出对象的数据成员.
*                                                                       *
*************************************************************************
*/

void DynamicsLevelIntegrator::printClassData(ostream& os) const
{
   os << "\nDynamicsLevelIntegrator::printClassData..." << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_cfl = " << d_cfl << endl;
   os << "d_initial_dt = " << d_initial_dt << endl;
}


