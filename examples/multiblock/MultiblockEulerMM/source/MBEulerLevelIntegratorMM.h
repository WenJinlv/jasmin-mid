 //
// 文件名:  MBEulerLevelIntegratorMM.h
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  移动网格方法求解Euler方程组的网格层时间积分算法
//

#ifndef included_MBEulerLevelIntegratorMM
#define included_MBEulerLevelIntegratorMM
 
#include "JASMIN_config.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "tbox/Pointer.h"

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "ReductionIntegratorComponent.h"

#include "MultiblockGridGeometry.h"

using namespace JASMIN;

#include "MBEulerMM.h"

/**
 * @brief 该类实现网格层时间积分算法策略类. 采用移动网格方法求解Euler方程.
 * 该类仅考虑单层网格.
 *
 * 该类实现如下抽象接口函数:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *    -# acceptTimeDependentSolution()
 *
 * 该类对象从输入文件或重启动数据库读取如下参数(输入文件优先):
 * - 选择读取的参数:
 *    - cfl \n
 *       双精度浮点型, 表示在时间步进的过程中, 控制时间步长的CFL条件. 缺省值为0.9.
 *
 * 该类的对象注册到重启动数据库.
 * 
 * @see algs::TimeIntegratorLevelStrategy, algs::IntegratorComponentManager
 *
 */

class MBEulerLevelIntegratorMM :
public tbox::Serializable,
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    */     
   MBEulerLevelIntegratorMM(const string& object_name,
           tbox::Pointer<tbox::Database> input_db,
           MBEulerMM* patch_strategy,
           tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
           const bool register_for_restart = true);

   /**
    * @brief 析构函数.
    */
   virtual ~MBEulerLevelIntegratorMM();

   /// @name 必须实现如下纯虚函数.
   //@{

   /**
    * @brief 初始化网格层时间积分算法. 创建所有积分构件.
    */
   void 
   initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 初始化指定网格层的数据片.
    * 初值构件类 algs::InitializeIntegratorComponent 支撑该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向待初始化网格层所在的网格片层次结构.
    * @param level_number 输入参数, 整型, 表示待初始化的网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示该网格层可被进一步细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param old_level 输入参数, 指针, 指向初始化网格层的旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示初始化为数据片调度内存空间. 
    *
    */
   void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                       const int    level_number,
                       const double init_data_time,
                       const bool   can_be_refined,
                       const bool   initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level =
                               tbox::Pointer< hier::BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /**
    * @brief 获取指定网格层的时间步长. 步长构件类 algs::DtIntegratorComponent
    * 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
    *
    * @return 双精度浮点型, 表示网格层的时间步长.
    *
    */
   double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层积分一个时间步. 数值构件类 algs::NumericalIntegratorComponent
    * 支撑该函数的实现.
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
    *
    * 断言检查: 指针level非空.
    */ 
   int
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
    * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
    * 
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    * 
    * 断言检查: 参数level非空.
    */
   void 
   acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                  const double new_time, 
                  const bool deallocate_data);

   // @}

   /**
    * @brief 输出数据成员到指定的数据流.
    *
    * @param os 输入参数, 输出流类, 用于存储输出的结果.
    */
   void printClassData(ostream& os) const;

   /**
    * @brief 输出数据成员到指定的数据库.
    *
    * @param db 输入参数, 指针, 指向存储输出结果的数据库.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);


private:

   /**
    * @brief 从输入数据库中读取参数, 初始化数据成员.
    *
    * @param db 输入参数, 指针, 指向输入数据库.
    * @param is_from_restart 输入参数, 逻辑型, 真值表示从重启动数据库读取参数.
    *
    * 断言检查: 参数db非空.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                             bool is_from_restart);

   /**
    * @brief 从重启动数据库读取参数, 初始化数据成员.
    */
   void getFromRestart();

   // 对象名称.
   string d_object_name;

   bool d_registered_for_restart; /*!< 真值表示注册该对象到重启动管理器. */
   double d_cfl;           /*!< Courant-Friedrichs-Levy常数. 缺省值为0.9. */

   // 网格片时间积分算法对象.
   MBEulerMM* d_patch_strategy;
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > d_grid_geometry;

   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >
                                                       d_init_intc;

   // 步长构件.
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> > d_dt_intc;

   // 数值构件: 计算梯度, 求Euler方程, 求控制函数, 移动网格.
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > 
                                                       d_grad_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                       d_pde_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                       d_wc0_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                       d_wci_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                       d_coords_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                       d_movegrad_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > 
                                                       d_movesolu_intc;

   // 复制构件: 接收数值解, 移动网格迭代.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >
                                                       d_reset_intc;
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >
                                                       d_cycle_copy_intc;

   // 内存构件: 新值, 演算. 
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_new_intc;
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_scratch_intc;

   tbox::Pointer< algs::ReductionIntegratorComponent<NDIM> >
                                                       d_iter_error_intc;
};

#endif
