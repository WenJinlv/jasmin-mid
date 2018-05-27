//
// 文件名:  DynamicLevelIntegrator.h
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  演示非协调拼接构件使用方法的网格层时间积分算法类.
//

#ifndef included_DynamicLevelIntegrator
#define included_DynamicLevelIntegrator

#include "JASMIN_config.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "tbox/Pointer.h"

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "IntegratorComponentManager.h"
#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "Dynamic.h"
#include "MultiblockGridGeometry.h"
#include "MultiblockPatchHierarchy.h"
#include "NonconformContactIntegratorComponent.h"

using namespace JASMIN;

/**
 * @brief 网格层时间积分算法类.
 *
 * 应用必须实现以下函数:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *    -# acceptTimeDependentSolution()
 *
 * 以上所有函数均可通过装配相应的积分构件实现.
 * 积分构件由构件管理器 algs::IntegratorComponentManager创建. 
 * 在各个函数的说明中, 将列出该函数可能使用的积分构件. 
 *
 * 该类对象从输入文件或重启动数据库读取如下参数(输入文件优先):
 * - 选择读取的参数:
 *    - cfl \n
 *       双精度浮点型, 表示在时间步进的过程中, 控制时间步长的CFL条件. 缺省值为0.9.
 *
 * 该类的对象被注册到重启动数据库.
 * 
 * @see algs::ExplicitTimeIntegrator, algs::IntegratorComponentManager
 *
 */

class DynamicLevelIntegrator :
      public virtual tbox::Serializable,
      public algs::TimeIntegratorLevelStrategy<NDIM>
  {
  public:
    /**
     * @brief 构造函数.
     */
    DynamicLevelIntegrator(const string& object_name,
                          tbox::Pointer<tbox::Database> input_db,
                          Dynamic* patch_strategy,
                          tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry,
                          tbox::Pointer< hier::MultiblockPatchHierarchy<NDIM> > patch_hierarchy,
                          const bool register_for_restart = true,
                          const bool use_time_refinement = true,
                          const bool use_implicit_discretization=false);

    /**
     * @brief 析构函数.
     */
    ~DynamicLevelIntegrator();

    /// @name 必须实现如下纯虚函数.
    //@{

    /**
     * @brief 初始化网格层时间积分算法. 
     *
     * 该函数根据实际应用和计算方法, 由积分构件管理器创建并初始化所有积分构件.
     */
    void
    initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

    /**
     * @brief 初始化指定网格层的数据片.
     *
     * @param hierarchy      输入参数, 指针,         指向待初始化网格层所在的网格片层次结构.
     * @param level_number   输入参数, 整型,         待初始化的网格层层号.
     * @param init_data_time 输入参数, 双精度浮点型, 初始化的时刻.
     * @param can_be_refined 输入参数, 逻辑型,       真值表示该网格层可被进一步细化.
     * @param initial_time   输入参数, 逻辑型,       真值表示当前时刻为计算的初始时刻.
     * @param old_level      输入参数, 指针,         指向初始化网格层的旧网格层.
     * @param allocate_data  输入参数, 逻辑型,       真值表示初始化为数据片调度内存空间. 
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
     * @param level        输入参数, 指针,         指向网格层.
     * @param dt_time      输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
     * @param initial_time 输入参数, 逻辑型,       真值表示当前时刻为计算的初始时刻.
     * @param flag_last_dt 输入参数, 整型,         表示上个时间步积分返回的状态.
     * @param last_dt      输入参数, 双精度浮点型, 表示上个时间步长.
     *
     * @return 双精度浮点型, 表示网格层的时间步长.
     */
    double
    getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
               const double dt_time,
               const bool initial_time,
               const int  flag_last_dt,
               const double last_dt);

    /**
     * @brief 在网格层上积分一个时间步.
     *
     * @param level                  输入参数, 指针,         指向待积分的网格层.
     * @param current_time           输入参数, 双精度浮点型, 当前时间步的起始时刻.
     * @param predict_dt             输入参数, 双精度浮点型, 为该时间步预测的时间步长.
     * @param max_dt                 输入参数, 双精度浮点型, 允许的最大时间步长.
     * @param min_dt                 输入参数, 双精度浮点型, 允许的最小时间步长.
     * @param first_step             输入参数, 逻辑型,       真值当前为重构后或时间步序列的第1步.
     * @param hierarchy_step_number, 输入参数, 整型,         表示网格片层次结构的积分步数,
     *                                                       也就是最粗网格层的积分步数.
     * @param actual_dt              输出参数, 双精度浮点型, 当前时间步进实际采用的时间步长.
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
     * @brief 在指定的网格层上, 时间步积分完毕, 存储最新的数值解, 
     * 更新网格层的状态到新的时刻, 完成下一个时间步积分的准备工作.
     * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
     * 
     * @param level           输入参数, 指针,         指向网格层.
     * @param new_time        输入参数, 双精度浮点型, 表示新的时刻.
     * @param deallocate_data 输入参数, 逻辑型,       真值表示接收数值解后, 释放新值数据片的内存空间.
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

    ///@name 私有成员函数.
    // @{

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

    //@}

    // 对象名称.
    string d_object_name;

    bool d_registered_for_restart; /*!< 真值表示注册该对象到重启动管理器. */

    // 网格片时间积分算法对象.
    Dynamic* d_patch_strategy;

    // 多块网格几何类对象.
    tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > d_grid_geometry;

    // 多块网格片层次结构.
    tbox::Pointer< hier::MultiblockPatchHierarchy<NDIM> > d_patch_hierarchy;

    // 初始化构件.
    tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >  d_init_component;

    // 步长构件.
    tbox::Pointer< algs::DtIntegratorComponent<NDIM> >          d_dt_component;

    // 数值构件: 计算下一时刻的网格结点位置坐标.
    tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_update_node_coordinate_component;

    // 复制构件: 将某些数据片从NONCONFORM上下文复制到CURRENT上下文.
    tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >     d_accept_component;

    // 内存构件: 管理NONCONFORM上下文中数据片的内存空间.
    tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_memory_component;

    // 非协调拼接构件: 验证非协调拼接构件的正确性.
    tbox::Pointer< algs::NonconformContactIntegratorComponent<NDIM> > d_nonconform_component;

    // 统计非协调拼接计算时间的计时器.
    tbox::Pointer<tbox::Timer> t_nonconform_contact_computing;
  };

#endif
