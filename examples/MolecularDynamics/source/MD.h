//
// 文件名: MD.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 (浜, 28  9 2007) $
// 描述  : 粒子模拟问题的网格片时间积分算法类的实现.
//

#ifndef included_MD
#define included_MD

#include "tbox/Tracer.h"
#include "JASMIN_config.h"
#include "IntVector.h"
#include "tbox/Pointer.h"
#include "ParticleVariable.h"
#include "IntegratorComponent.h"
#include "tbox/Serializable.h"
#include "tbox/Database.h"
#include "tbox/Array.h"

#include "Patch.h"
#include "CellVariable.h"
#include "FaceVariable.h"
#include "VariableContext.h"

#include "CartesianCoordinates.h"
#include "UniRectangularGridGeometry.h"
#include "StandardComponentPatchStrategy.h"

#include "JaVisDataWriter.h"

using namespace JASMIN;

/**
 * @brief 围绕分子动力学问题,
 * 实现网格片时间积分算法基类 algs::StandardComponentPatchStrategy.
 *
 * 必须实现如下纯虚函数:
 *  - initializeComponent() : 初始化积分构件.
 *  - computeOnPatch() : 为数值构件完成数值计算.
 *  - getPatchDt()     :  为步长构件计算时间步长.
 *  - initializePatchData() : 为初值构件完成初始化.
 *
 * 数据通信时, 如果需要为数据片填充物理边界条件, 则必须重载如下函数:
 *  - setPhysicalBoundaryConditions() : 填充物理边界条件.
 *
 * @see algs::IntegratorComponent
 */

class MD :
            public tbox::Serializable,
            public algs::StandardComponentPatchStrategy<NDIM>
{
public:


    /*! @brief 构造函数.
     * @param object_name 输入参数, 字符串, 表示对象名称.
     * @param input_db    输入参数, 指针, 指向输入数据库.
     * @param grid_geom   输入参数, 指针, 指向均匀矩形网格几何.
     * @note
     * 该函数执行以下操作:
     *  -# 初始化数据成员;
     *  -# 创建定义物理模型控制方程中的数值解的变量;
     *  -# 如果需断点续算, 则从重启动文件中读入数据;
     *  -# 从输入数据库中读入数据, 该操作可能覆盖从重启动文件中读入的数据.
     *  -# 将当前对象注册到重启动管理器.
     *
     */
    MD(const string& object_name,
       tbox::Pointer<tbox::Database> input_db,
       tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom);

    /*! @brief构函数.
     *
     * @note 虚的空函数.
     */
    virtual ~MD();

    /*!
     * @brief 初始化指定的积分构件, 
     *        注册待填充的数据片或待调度内存空间的数据片到积分构件.
     *
     * @param intc 输入参数, 指针, 指向待初始化的积分构件对象.
     *
     * 断言检查指针非空.
     */
    virtual void initializeComponent(
        algs::IntegratorComponent<NDIM>* intc) const;


    /*!
     * @brief 当前积分构件为数值构件. 该函数完成数值计算.
     *
     * @param patch 输入参数, 网格片类, 表示网格片.
     * @param time  输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt    输入参数, 双精度浮点型, 表示时间步长.
     * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
     * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
     */
    virtual void computeOnPatch(hier::Patch<NDIM>& patch,
                                const double  time,
                                const double  dt,
                                const bool    initial_time,
                                const string& intc_name);

    /**
     * @brief 当前积分构件为初值构件. 该函数初始化数据片.
     *
     * @param patch 输入参数, 网格片类, 表示网格片.
     * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
     * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
     * @param intc_name 输入参数, 字符串, 表示初值构件的名称.
     */
    virtual void initializePatchData(
        hier::Patch<NDIM>& patch,
        const double  time,
        const bool    initial_time,
        const string& intc_name);

    /**
     * @brief 计算粒子受力，更新粒子位置.
     *
     * @param patch 输入参数, 网格片类, 表示网格片.
     * @param current_time  输入参数, 双精度浮点型, 表示计算时刻.
     * @param dt            输入参数, 双精度浮点型, 表示时间步长.
     */
    void updateParticleOnPatch(
        hier::Patch<NDIM>& patch,
        const double current_time,
        const double dt);

    /**
     * @brief 计算网格片单元上的负载.
     *
     * @param patch          输入参数, 网格片类, 表示网格片.
     * @param measured_load  输入参数, 双精度型, 表示本地处理器的当前负载.
     * @param numeric_step   输入参数, 整型, 表示当前负载调整包含的当前时间步数.
     * @param workload_index 输入参数, 整型, 表示存储负载的数据片索引号. 
     */
    void computeWorkloadOnPatch(
        hier::Patch<NDIM>& patch,
        const double measured_load,
        const int    numeric_step,
        const int    workload_index);

    /**
     * @brief 设置当前的时间步数.
     * @param integrator_step   输入参数, 整型, 表示当前时间步数.
     */
    void setIntegratorStep(const int integrator_step);

    /**
     * @brief 在积分构件中, 为数据片填充物理边界条件.
     *
     * @param patch     输入参数, 网格片类, 表示网格片.
     * @param fill_time 输入参数, 双精度浮点型, 表示当前时刻.
     * @param ghost_width_to_fill 输入参数, 整型向量, 表示物理边界影像区被填充的宽度.
     * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
     */
    virtual void setPhysicalBoundaryConditions(
        hier::Patch<NDIM>& patch,
        const double fill_time,
        const hier::IntVector<NDIM>& ghost_width_to_fill,
        const string& intc_name);

   /*! 
    * @brief 当前积分构件为步长构件. 该函数计算时间步长.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt   输入参数, 整型, 表示前次积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示前次积分的时间步长.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    *
    * @return 双精度浮点型, 表示时间步长.
    */

    virtual double getPatchDt(hier::Patch<NDIM>& patch,
                              const double  time,
                              const bool    initial_time,
                              const int     flag_last_dt,
                              const double  last_dt,
                              const string& intc_name);


    /*!
     * @brief 将数据成员输出到重启动数据库.
     * @param db 输入参数, 指针, 指向重启动数据库.
     * @note
     * 该函数是对类 JASMIN::tbox::Serializable 中定义的函数接口的实现.
     */
    void putToDatabase(tbox::Pointer<tbox::Database> db);

    /*!
     * @brief 注册JaVis数据输出器到JASMIN框架, 同时, 注册可视化数据片.
     * @param viz_writer 输入参数, 指针, 表示 JaVis 数据输出器.
     */
    void registerJaVisDataWriter(
        tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer);

    /*!
     * @brief 输出当前对象中所有数据成员到指定输出流.
     * @param os 输入输出参数, 输出流类型, 表示指定的输出流.
     */
    void printClassData(ostream& os) const;

    /*!
     * @brief 当前积分构件为归约构件. 该函数计算粒子总数目.
     *
     * @param patch 输入参数, 网格片类, 表示网格片.
     * @param time  输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt    输入参数, 双精度浮点型, 表示时间步长.
     * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
     * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
     */
    virtual void reduceOnPatch(double *vector,
                             int     len,
                             hier::Patch<NDIM> &patch,
                             const double time,
                             const double dt,
                             const string& intc_name);

private:

    /*
     * 从输入数据库和重启动数据库读入数据, 但是, 输入数据库优先于从启动数据库.
     * 断言检查不允许数据库指针为空.
     */
    void getFromInput(tbox::Pointer<tbox::Database> db,
                      bool is_from_restart);

    void getFromRestart();

    // 创建所有变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
    void registerModelVariables();


    /*
     * 对象名.
     */
    string d_object_name;

    /*
     * 坐标系和网格几何.
     */
    tbox::Pointer<geom::CartesianCoordinates<NDIM> >   d_coords;
    tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > d_grid_geometry;

    /*
     * 可视化数据输出器.
     */
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > d_Javis_writer;

    /**
     * 存储数值解 - [u]
     */
    tbox::Pointer< pdat::ParticleVariable<NDIM,double> > d_Cu_particle_var;
    tbox::Pointer< pdat::CellVariable<NDIM, double> > d_particle_density_var;
    tbox::Pointer< pdat::CellVariable<NDIM, double> > d_patch_mapping_var;

    // 上下文: 当前值, 新值, 可视化.
    tbox::Pointer< hier::VariableContext > d_current;
    tbox::Pointer< hier::VariableContext > d_new;
    tbox::Pointer< hier::VariableContext > d_plot_context;

    // 数据片索引号.
    int d_Cu_particle_id;
    int d_particle_density_id;
    int d_patch_mapping_id;

    int d_average_number_in_cell;
    double d_particle_memory_grow;

    hier::IntVector<NDIM> d_nghosts;

    string d_model_problem;
    int d_integrator_step;
    /*
    * 说明问题类型和初始条件的变量.
    */
    string d_data_problem;
    int d_data_problem_int;
    double d_time_step;
    int d_ipt_cell;
    double d_triangle_vertex;
    double d_triangle_angle;
    double d_triangle_depth;

    // 物理模型参数. 左边物体的速度和相对位置, 右边物体的速度和相对位置.
    double  d_left_obj_velocity, d_right_obj_velocity;
    double d_left_obj_position[2], d_right_obj_position[2];

    /*
     * 统计量.
    */
    //tbox::Pointer<tbox::Statistic> d_number_particles_on_patch;

};

#endif
