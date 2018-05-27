//
// 文件名: Dynamic.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date$
// 描述  : 演示非协调拼接构件使用方法的网格片时间积分算法类.
//

#ifndef included_Dynamic
#define included_Dynamic

#include "JASMIN_config.h"
#include "Index.h"
#include "IntVector.h"
#include "tbox/Pointer.h"
#include "IntegratorComponent.h"
#include "tbox/Serializable.h"
#include "tbox/Database.h"
#include "tbox/Array.h"

#include "Patch.h"
#include "CellVariable.h"
#include "NodeVariable.h"
#include "SideVariable.h"
#include "VariableContext.h"

#include "CartesianCoordinates.h"
#include "MultiblockDeformingGridGeometry.h"
#include "StandardComponentPatchStrategy.h"

#include "JaVisDerivedDataStrategy.h"
#include "JaVisDataWriter.h"

#include "DeformingGridInputUtilities2.h"

// 表示重启动文件的版本号的常量
#define Dynamic_version (1)

using namespace JASMIN;

/**
 * @brief 类 Dynamic 提供函数接口, 实现初始化, 计算时间步长以及各数值模块的计算功能. 
 *        各函数调用Fortran子程序完成计算.
 * 
 * 必须实现如下纯虚函数:
 *  - initializeComponent(): 初始化积分构件.
 *  - initializePatchData(): 为初始化构件完成初始化操作.
 *  - getPatchDt():          为步长构件计算时间步长.
 *  - computeOnPatch():      为数值构件完成数值计算.
 *
 * 细化数据通信时, 如果需要为数据片填充物理边界条件, 则必须重载如下函数:
 *  - setPhysicalBoundaryConditions(): 填充物理边界条件.
 *
 * @see algs::IntegratorComponent
 */

class Dynamic :
      public tbox::Serializable,
      public algs::StandardComponentPatchStrategy<NDIM>
  {
  public:

    /*! @brief 构造函数.
     * @param object_name 输入参数, 字符串, 表示对象名称.
     * @param input_db    输入参数, 指针,   指向输入数据库.
     * @param grid_geom   输入参数, 指针,   指向均匀矩形网格几何.
     * @note
     * 该函数执行以下操作:
     *  -# 初始化数据成员;
     *  -# 创建变量;
     *  -# 如果需断点续算, 则从重启动文件中读入数据;
     *  -# 从输入数据库中读入数据, 该操作可能覆盖从重启动文件中读入的数据.
     *  -# 将当前对象注册到重启动管理器.
     *
     */
    Dynamic(const string& object_name,
            tbox::Pointer<tbox::Database> input_db,
            tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom,
            tbox::Pointer<appu::DeformingGridInputUtilities2> grid_tool,
            int module_id);

    /*! @brief 析构函数.
     *
     * @note 虚的空函数.
     */
    ~Dynamic();

    /*!
     * @brief 初始化指定的积分构件, 
     *        注册待填充的数据片或待调度内存空间的数据片到积分构件.
     *
     * @param component 指针, 输入参数, 指向待初始化的积分构件对象.
     *
     * 断言检查指针非空.
     */
    void initializeComponent(
      algs::IntegratorComponent<NDIM>* component) const;

    /*!
     * @brief 该函数在单个网格片上计算时间步, 服务于步长构件. 
     *
     * @param patch          输入参数, 网格片类,     表示网格片.
     * @param time           输入参数, 双精度浮点型, 当前时刻.
     * @param initial_time   输入参数, 逻辑型,       真值表示当前为计算的初始时刻.
     * @param flag_last_dt   输入参数, 整型,         前次积分返回的状态.
     * @param last_dt        输入参数, 双精度浮点型, 前次积分的时间步长.
     * @param component_name 输入参数, 字符串,       积分构件的名称.
     *
     * @return 双精度浮点型, 表示时间步长.
     */
    double getPatchDt(hier::Patch<NDIM>& patch,
                              const double  time,
                              const bool    initial_time,
                              const int     flag_last_dt,
                              const double  last_dt,
                              const string& component_name);

    /*!
     * @brief 该函数完成数值计算, 服务于数值构件.
     *
     * @param patch          输入参数, 网格片类,     表示网格片.
     * @param time           输入参数, 双精度浮点型, 当前时刻.
     * @param dt             输入参数, 双精度浮点型, 时间步长.
     * @param initial_time   输入参数, 逻辑型,       真值表示当前为计算的初始时刻.
     * @param component_name 输入参数, 字符串,       积分构件的名称.
     */
    void computeOnPatch(hier::Patch<NDIM>& patch,
                                const double  time,
                                const double  dt,
                                const bool    initial_time,
                                const string& component_name);

    /*!
     * @brief 该函数完成非协调拼接计算.
     *
     * @param patch          输入参数, 网格片类,     表示网格片.
     * @param location       输入参数, 整型,         非协调块边界与当前网格片的相对位置, 详细定义见注释.
     * @param ghost_patch    输入参数, 网格片类,     与当前网格片相应的非协调影像网格片.
     * @param ghost_location 输入参数, 整型,         非协调块边界与非协调影像网格片的相对位置, 详细定义见注释.
     * @param time           输入参数, 双精度浮点型, 当前时刻.
     * @param dt             输入参数, 双精度浮点型, 时间步长.
     * @param component_name 输入参数, 字符串,       表示积分构件的名称.
     *
     * @note
     *  参数(ghost_)location的取值定义了非协调块边界与网格片的相对位置, 具体如下:
     *  - 0: 非协调块边界在网格片第0维方向的左边.
     *  - 1: 非协调块边界在网格片第0维方向的右边.
     *  - 2: 非协调块边界在网格片第1维方向的左边.
     *  - 3: 非协调块边界在网格片第1维方向的右边.
     *  - 4: 非协调块边界在网格片第2维方向的左边.
     *  - 5: 非协调块边界在网格片第2维方向的右边.
     */
    void computeNonconformContactOnPatch(hier::Patch<NDIM>& patch,
        const int location,
        hier::Patch<NDIM>& ghost_patch,
        const int ghost_location,
        const int bdry_type,
        const double time,
        const double dt,
        const string& component_name);

    /**
     * @brief 该函数初始化数据片, 服务于初始化构件.
     *
     * @param patch          输入参数, 网格片类,     表示网格片.
     * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
     * @param initial_time   输入参数, 逻辑型,       真值表示当前为计算的初始时刻.
     * @param component_name 输入参数, 字符串,       初始化构件的名称.
     */
    void initializePatchData(hier::Patch<NDIM>& patch,
                                     const double  time,
                                     const bool    initial_time,
                                     const string& component_name);


    /**
    * @brief 设置物理边界条件.
    *
    * @param patch               输入参数, 网格片类,     表示网格片.
    * @param fill_time           输入参数, 双精度浮点型, 表示当前时刻.
    * @param ghost_width_to_fill 输入参数, 整型向量,     表示物理边界影像区被填充的宽度.
    * @param component_name      输入参数, 字符串,       积分构件的名称.
    */
    void setPhysicalBoundaryConditions(
      hier::Patch<NDIM>& patch,
      const double fill_time,
      const hier::IntVector<NDIM>& ghost_width_to_fill,
      const string& component_name);

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
    void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer);

    /*!
     * @brief 输出当前对象中所有数据成员到指定输出流.
     * @param os 输入输出参数, 输出流类型, 表示指定的输出流.
     */
    void printClassData(ostream& os) const;

    /**
     * 将参数patch指定网格片上的数据片以参数prec指定的精度输出到参数os指定的
     * 输出流中. 其中, 实参prec为数据的打印精度(十进制位数), -1表示使用缺省
     * 精度(例如, double类型和complex类型为12位, float类型为6位).
     */

    void printPatchData(hier::Patch<NDIM>& patch,
                        const int patchdata_id,
                        const int patchdata_type,
                        const int patchdata_depth,
                        ostream& os = tbox::plog,
                        int prec = -1);

  private:

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
     * 模型编号.
     */
    int d_module_id;

    /*
     * 网格几何.
     */
    tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > d_grid_geometry;

    /*
     * 二维变形网格辅助生成工具.
     */
    tbox::Pointer< appu::DeformingGridInputUtilities2 > d_grid_tool;

    /*
     * 可视化数据输出器.
     */
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > d_javis_writer;

    // 速度, 定义在网格结点上, 深度为2(2d), 3(3d).
    tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_velocity;
    // 网格结点坐标.
    tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_node_coord;

    // 数据片索引号.
    // 速度.
    int d_velocity_current_id;
    int d_velocity_nonconform_id;
    // 网格结点坐标.
    int d_node_coord_current_id;
    int d_node_coord_nonconform_id;

    // 非协调块边界数组.
    tbox::Array<int> d_nonconform_boundaries;

    // 控制参数.
    double d_small;
    double d_large;
    // 时间步长.
    double d_dt;

    // 网格块移动速度.
    tbox::Array<double> d_velocity_in;

    // 处理非协调块边界的计算格式的宽度.
    int d_opstencil_cells;

    // 警戒单元数目.
    int d_guard_cells;

    // 上下文: 当前值(current), 非协调接触(nonconform), 可视化(plot).
    tbox::Pointer< hier::VariableContext > d_current_context;
    tbox::Pointer< hier::VariableContext > d_nonconform_context;
    tbox::Pointer< hier::VariableContext > d_plot_context;
  };

#endif
