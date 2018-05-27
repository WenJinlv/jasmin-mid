 //
// 文件名:  EulerLevelIntegrator.h
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  Euler方程组的网格层时间积分算法
//

#ifndef included_EulerLevelIntegrator
#define included_EulerLevelIntegrator
 
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

using namespace JASMIN;

/**
 * @brief 该类以显式时间离散格式实现网格层时间积分算法策略类.
 *
 * 单层网格的应用必须实现:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *    -# acceptTimeDependentSolution()
 *
 * 多层网格的应用必须实现:
 *    -# usingRefinedTimestepping()
 *    -# synchronizeNewHierarchy()
 *    -# synchronizeCoarserLevel()
 *    -# tagCellsForRefinement()
 *
 * 以上所有函数均使用了积分构件. 在各个函数的说明中,
 * 将列出该函数使用的积分构件的类型. 
 *
 * 该类对象从输入文件或重启动数据库读取如下参数(输入文件优先):
 * - 选择读取的参数:
 *    - cfl \n
 *       双精度浮点型, 表示在时间步进的过程中, 控制时间步长的CFL条件. 缺省值为0.9.
 *
 * 该类的对象注册到重启动数据库.
 * 
 * @see algs::TimeIntegratorLevelStrategy
 *
 */

class EulerLevelIntegrator :
public tbox::Serializable,
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    */     
   EulerLevelIntegrator(const string& object_name,
                   tbox::Pointer<tbox::Database> input_db,
                   algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
                   tbox::Pointer< hier::GridGeometry<NDIM> > grid_geometry,
                   const bool register_for_restart = true,
                   const bool use_time_refinement = true);

   /**
    * @brief 析构函数.
    */
   virtual ~EulerLevelIntegrator();

   /// @name 必须实现如下纯虚函数.
   //@{

   /**
    * @brief 初始化网格层时间积分算法. 
    *
    * 该函数根据实际应用和计算方法, 创建所有积分构件.
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

   /// @name 如果当前计算为显式时间积分且网格片层次结构包含多层网格, 则须重载如下虚函数:
   // @{

   /**
    * @brief 网格片层次结构被创建或重建后, 从细网格层到粗网格层, 依次同步.
    * 同步构件 algs::SynchronizeIntegratorComponent 支持该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向新创建或重建的网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示待同步的最粗网格层层号.
    * @param finest_level 输入参数, 整型, 表示待同步的最细网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前为计算的初始时刻.
    */
   virtual void synchronizeNewHierarchy(
                   const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                   const int coarsest_level, 
                   const int finest_level, 
                   const double sync_time, 
                   const bool initial_time);

   /**
    * @brief 时间积分的步进模式.
    * @return 真值表示细化模式, 假值表示同步模式.
    */
   virtual bool
   usingRefinedTimestepping();

   /**
    * @brief 细网格层同步相邻的粗网格层.
    * 同步构件 algs::SynchronizeIntegratorComponent 支持该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param finer_level_number 输入参数, 双精度型, 表示细网格层的层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步的时刻.
    * @param coarser_old_time  输入参数, 双精度型, 表示粗网格层当前时间步的起始时刻.
    */
   virtual void 
   synchronizeCoarserLevel(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int finer_level_number,
      const double sync_time, 
      const double coarser_old_time);

   /*!
    * @brief 在指定的网格层上, 标记待细化的网格单元. 
    * 数值构件 algs::NumericalIntegratorComponent 支持该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param level_number 输入参数, 整型, 网格层的层号.
    * @param error_data_time 输入参数, 双精度浮点型, 估计误差的时刻.
    * @param tag_index  输入参数, 整型, 存储标记值的网格数据片的索引号.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为初始时刻.
    *
    * @note
    *  目前, 参数regrid_start_time仅在Richard外插算中使用.
    *  如果误差粗化率为2, 则该参数值为前一时间步的起始时刻; 如果误差粗化率为3,
    *  则该参数值为前两个时间步的起始时刻.
    */
   virtual void
   tagCellsForRefinement(const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
                         const int level_number,
                         const double error_data_time,
                         const int tag_index,
                         const bool initial_time);

   //@}

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


   ///@name 私有成员函数.
   // @{
      
      /**
       *  @brief 标准时间积分函数, 在时间步序列中, 积分1个时间步.
       *  @param level        输入参数, 指针, 指向网格层.
       *  @param current_time 输入参数, 双精度浮点型, 表示积分的起始时刻.
       *  @param actual_dt    输入参数, 双精度浮点型, 表示积分步长.
       *  @param first_step   输入参数, 逻辑型, 真值表示当前为时间步序列的第1个时间步.
       */
      void standardAdvanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                const double current_time,
                const double actual_dt,
                const bool   first_step);

      /**
       *  @brief 同步时间积分函数, 为了配合细网格层对粗网格层的校正, 积分1个时间步.
       *  @param coarser_level 输入参数, 指针, 指向粗网格层.
       *  @param finer_level   输入参数, 指针, 指向细网格层.
       *  @param current_time  输入参数, 双精度浮点型, 表示积分的起始时刻.
       *  @param actual_dt     输入参数, 双精度浮点型, 表示积分步长.
       */
      void synchronizeAdvanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > coarser_level,
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > finer_level,
                const double current_time,
                const double actual_dt);

   //@}

   // 对象名称.
   string d_object_name;

   bool d_use_time_refinement;    /*!< 真值表示采用时间细化模式. */
   bool d_registered_for_restart; /*!< 真值表示注册该对象到重启动管理器. */
   double d_cfl;                  /*!< Courant-Friedrichs-Levy常数. 缺省值为0.9. */

   // 网格片时间积分算法对象.
   algs::StandardComponentPatchStrategy<NDIM>* d_patch_strategy;
   tbox::Pointer< hier::GridGeometry<NDIM> > d_grid_geometry;

   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >  d_init_intc;

   // 步长构件.
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> >          d_dt_intc;

   // 数值构件:计算通量, 守恒差分, 标记待细化的网格单元.
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >   d_flux_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >   d_diff_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >   d_tag_intc;

   // 同步构件: 粗细网格层交接面流连续, 细网格层覆盖的粗网格层内部区域.
   tbox::Pointer< algs::SynchronizeIntegratorComponent<NDIM> > d_sync_flux_intc;
   tbox::Pointer< algs::SynchronizeIntegratorComponent<NDIM> > d_sync_overlay_intc;

   // 复制构件: 接收数值解, 同步复制.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >     d_reset_intc;
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >     d_sync_copy_intc;

   // 内存构件: 新值数据片, 演算数据片. 
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_new_intc;
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_scratch_intc;

};

#endif
