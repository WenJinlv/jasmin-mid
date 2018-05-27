//
// 文件名: MDLevelIntegrator.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 170 $
// 修改  : $Date: 2007-06-27 08:44:22 +0800 (涓, 27  6 2007) $
// 描述  : 粒子模拟问题的网格层时间积分算法的实现.
//

#ifndef included_MDLevelIntegrator
#define included_MDLevelIntegrator

#include "JASMIN_config.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Tracer.h"

#include "ParticleVariable.h"
#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "ParticleCommComponent.h"
#include "ReductionIntegratorComponent.h"
#include "MD.h"

using namespace JASMIN;


/**
 * @brief 该类以显式时间离散格式实现网格层时间积分算法策略类. 
 * 
 * 单层网格的应用必须实现:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *
 * 该类的对象注册到重启动数据库.
 * 
 * @see algs::ExplicitTimeIntegrator
 *
 */

class MDLevelIntegrator :
      public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
    /**
     * @brief 构造函数.
     *
     * @param object_name 输入参数, 字符串, 表示对象的名字.
     * @param input_db    输入参数, 指针, 指向输入数据库.
     * @param patch_strategy 输入参数, 指针, 指向网格片时间积分策略类.
     * @param md_model       输入参数, 指针, 指向分子动力学模型类.
     * @param grid_geometry  输入参数, 指针, 指向网格几何.
     *
     * 断言检查: object_name, input_db, patch_strategy非空.
     */  
    MDLevelIntegrator(const string& object_name,
                      tbox::Pointer<tbox::Database> input_db,
                      algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
                      MD* md_model,
                      tbox::Pointer< hier::GridGeometry<NDIM> > grid_geometry);

    /**
     * @brief 析构函数.
     */
    virtual ~MDLevelIntegrator();

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
     * @brief 计算网格片中单元上的负载.
     *
     * @param patch_level 输入参数, 指针, 指向网格层.
     * @param measured_load  输入参数, 双精度型数组, 表示本地处理器的当前负载.
     * @param numeric_step   输入参数, 整型, 表示当前负载调整包含的当前时间步数.
     * @param workload_index 输入参数, 整型, 表示存储负载的数据片索引号. 
     */

    virtual void refreshNonUniformWorkload(
        const tbox::Pointer< hier::BasePatchLevel<NDIM> > patch_level,
        const tbox::Array<double>& measured_load,
        const int    numeric_step,
        const int    workload_index);

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
     * @brief 输出数据成员到指定的数据流.
     * @param os 输入参数, 数据流类型, 用于存储输出的结果.
     */
    virtual void printClassData(ostream& os) const;

    /**
     * 空函数, 不执行任何操作.
     */
    virtual void
    acceptTimeDependentSolution(
        const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
        const double new_time,
        const bool deallocate_data);

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

    void standardAdvanceLevel(
        const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
        const double current_time,
        const double actual_dt,
        const bool   first_step);


    /**
     * @brief 从重启动数据库读取参数, 初始化数据成员.
     */
    virtual void getFromRestart();

    string d_object_name;                       /*!< 对象名称. */

    double d_time_step;              /*!< 时间步长. */
    algs::StandardComponentPatchStrategy<NDIM>* d_patch_strategy;
    tbox::Pointer< hier::GridGeometry<NDIM> > d_grid_geometry;

    MD* d_md_model;

    // 初值构件.
    tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >  d_init_intc;
    // 粒子量通信构件: 迁移粒子、填充粒子 
    tbox::Pointer< algs::ParticleCommComponent<NDIM> >   d_particle_comm_intc;
    // 数值构件: 更新粒子位置.
    tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >   d_particle_comp_intc;
    // 规约构件：计算粒子总数目.
    tbox::Pointer< algs::ReductionIntegratorComponent<NDIM> >   d_particle_reduce_intc;

    tbox::Pointer< algs::DtIntegratorComponent<NDIM> >     d_dt_intc;

};

#endif


