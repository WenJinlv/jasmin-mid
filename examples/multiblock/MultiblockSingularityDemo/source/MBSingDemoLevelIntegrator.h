 //
// 文件名:  MBSingDemoLevelIntegrator.h
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  网格层时间积分算法.
//

#ifndef included_MBSingDemoLevelIntegrator
#define included_MBSingDemoLevelIntegrator
 
#include "JASMIN_config.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "tbox/Pointer.h"

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "OuterdataOperationIntegratorComponent.h"

#include "MultiblockGridGeometry.h"

using namespace JASMIN;

/**
 * @brief 该类实现网格层时间积分算法策略类.
 *
 * 该类实现如下抽象接口函数:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *    -# acceptTimeDependentSolution()
 *
 * @see algs::TimeIntegratorLevelStrategy
 */

class MBSingDemoLevelIntegrator :
public tbox::Serializable,
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    */     
   MBSingDemoLevelIntegrator(const string& object_name,
           tbox::Pointer<tbox::Database> input_db,
           algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
           tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
           const bool register_for_restart = true);

   /**
    * @brief 析构函数.
    */
   virtual ~MBSingDemoLevelIntegrator();

   /// @name 必须实现如下纯虚函数.
   //@{

   /**
    * @brief 初始化网格层时间积分算法. 创建所有积分构件.
    */
   virtual void 
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
   virtual void
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
   virtual double
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
    * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
    * 
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    * 
    * 断言检查: 参数level非空.
    */
   virtual void 
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
   virtual void printClassData(ostream& os) const;

   /**
    * @brief 输出数据成员到指定的数据库.
    *
    * @param db 输入参数, 指针, 指向存储输出结果的数据库.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);


private:

   /**
    * @brief 从输入数据库中读取参数, 初始化数据成员.
    *
    * @param db 输入参数, 指针, 指向输入数据库.
    * @param is_from_restart 输入参数, 逻辑型, 真值表示从重启动数据库读取参数.
    *
    * 断言检查: 参数db非空.
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db,
                             bool is_from_restart);

   /**
    * @brief 从重启动数据库读取参数, 初始化数据成员.
    */
   virtual void getFromRestart();

   // 对象名称.
   string d_object_name;

   bool d_registered_for_restart; /*!< 真值表示注册该对象到重启动管理器. */

   // 网格片时间积分算法对象.
   algs::StandardComponentPatchStrategy<NDIM>* d_patch_strategy;
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > d_grid_geometry;

   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >
                                                       d_init_intc;

   // 步长构件.
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> > d_dt_intc;

   // 数值构件.
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_sum_intc;

   // 复制构件.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> > d_reset_intc;

   // 内存构件: 新值, 演算. 
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_new_intc;
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_scratch_intc;

   // 外表面结点求累加和构件.
   tbox::Pointer<algs::OuterdataOperationIntegratorComponent<NDIM> > d_onsum_intc;

};

#endif
