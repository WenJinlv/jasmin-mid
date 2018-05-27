//
// 文件  :	DynamicsLevelIntegrator.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 1.12 $
// 修改  :	$Date: 2006/08/07 03:15:54 $
// 描述  :	程序流体力学在单层网格上的时间积分算法
//

#ifndef included_DynamicsLevelIntegrator
#define included_DynamicsLevelIntegrator

#ifndef included_JASMIN_config
#include "JASMIN_config.h"
#endif
#ifndef included_algs_TimeIntegratorLevelStrategy
#include "TimeIntegratorLevelStrategy.h"
#endif
#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "OuterdataOperationIntegratorComponent.h"
using namespace JASMIN;


#ifndef included_Dynamics
#include "Dynamics.h"
#endif

/**
 * @brief 流体力学网格层时间积分算法类 DynamicsLevelIntegrator
 * 是网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 的派生类.
 *
 * 该类通过调用用户网格片策略类 Dynamics, 
 * 以两层显式格式,将质量、动量和能量守恒方程从当前时刻t^n=current_time,
 * 积分一个Lagrange时间步, 到新的时刻t^{n+1}=new_time.
 * 当前,该类仅支持二维曲线自适应结构网格的单层网格.
 *
 * 该类需要从输入数据库读取如下参数. 
 *
 *    - \b cfl \n
 *         双精度浮点数,代表显式积分时间步长的Courant-Friedrichs-Levy条件.
 *         省缺值为1.0.      
 * 
 *    - \b initial_dt \n
 *         双精度浮点数,初始时间步长.
 */

class DynamicsLevelIntegrator :
public tbox::Serializable,
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    *
    * 从输入数据库或者重启动数据库初始化数据成员.
    */  
   DynamicsLevelIntegrator(
      const string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      Dynamics* patch_strategy,
      tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,  
      bool register_for_restart = true);

   /**
    * @brief 析构函数.
    */
   virtual ~DynamicsLevelIntegrator();

   /// @name 重载基类 algs::TimeIntegratorLevelStrategy 的成员函数:
   //@{

   /**
    * @brief 初始化网格层时间积分算法. 创建所有积分构件.
    */
   virtual void 
   initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 初始化网格层数据片.
    *
    * @param hierarchy 输入参数, 指针, 指向待初始化网格层所在的网格片层次结构.
    * @param level_number 输入参数, 整型, 表示待初始化的网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示该网格层可被进一步细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param old_level 输入参数, 指针, 指向初始化网格层的旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示初始化为数据片调度内存空间.
    */
   virtual void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
                       const int    level_number,
                       const double init_data_time,
                       const bool   can_be_refined,
                       const bool   initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level, 
                       const bool allocate_data );

   /**
    * @brief 获取网格层时间步长. 
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
    * @return 双精度浮点型, 表示网格层的时间步长.
    */
   virtual double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
              const double dt_time,
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层积分一个时间步. 
    *
    * 数值构件类 algs::NumericalIntegratorComponent 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向待积分的网格层.
    * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
    * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
    * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
    * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
    * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
    * @param hierarchy_step_number, 输入参数, 整型, 表示网格片层次结构的积分步数,
    *                               也就是最粗网格层的积分步数.
    * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
    *
    * @return 整型, 表示该时间步积分的状态, 其值的含义由用户负责解释.
    */
   virtual int
   advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt);
   
   /**
    * @brief 在指定的网格层上, 时间步积分完毕, 存储最新的获得数值解,
    * 更新网格层的状态到新的时刻, 完成下一个时间步积分的准备工作.
    *
    * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    */
   virtual void
   acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                  const double new_time,
                  const bool deallocate_data);

   //@}

   /** 
    * 将当前对象的数据成员输出到参数db指定的数据库.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

   /**
    * 将当前对象的数据成员输出到指定的输出流.
    */
   virtual void printClassData(ostream& os) const;


private:
   /**
    * 从指定的数据库读取参数初始化当前对象的数据成员.
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db,
                             bool is_from_restart);

   /**
    * 从重启动数据库读取数据初始化当前对象的数据成员. 
    */
   virtual void getFromRestart();


   /**
    * 流体力学用户网格片策略类对象,为本类提供操作数据片的数值函数.
    */
   Dynamics* d_patch_strategy;

   /**
    * 对象名
    */
   string d_object_name;

   /**
    * 逻辑变量, 用于控制是否向重启动数据库写入变量.
    */
   bool d_registered_for_restart;

   /**
    * 计算时间步长的控制参数.
    */
   double d_cfl;                  // Courant-Friedrichs-Levy限制条件
   double d_initial_dt;                  // first time step 


   // 多块网格几何.
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > d_grid_geometry;

   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >
                                                       d_init_intc;

   // 步长构件.
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> > d_dt_intc;

   // 数值构件: 前处理
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_preprocess_intc;
 
   // 数值构件：计算
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_velocity_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_grid_dens_energy_intc;

   // 数值构件：计算边压力
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_side_pressure_intc;

   // 外表面结点舍入误差构件.
   tbox::Pointer< algs::OuterdataOperationIntegratorComponent<NDIM> > 
                                                       d_outernode_roundoff;
   // 复制构件.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> > d_reset_intc;

   // 内存构件: 新值, 演算.
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_new_intc;
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_scratch_intc;


};

#endif
