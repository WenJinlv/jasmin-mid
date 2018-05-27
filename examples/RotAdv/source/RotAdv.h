//
// 文件名: RotAdv.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 (浜, 28  9 2007) $  
// 描述  : 刚体旋转问题的网格片时间积分算法类的实现.
//

#ifndef included_RotAdv
#define included_RotAdv

#include "JASMIN_config.h"
#include "IntVector.h"
#include "tbox/Pointer.h"
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
 * @brief 围绕刚体旋转问题,
 * 实现网格片时间积分算法基类 algs::StandardComponentPatchStrategy.
 *
 * 在均匀矩形的单网格片上, 求解刚体旋转对流方程:\n
 * @f$
 *     \frac{\partial u}{\partial t} + div(a*u) = 0
 * @f$
 * 其中, \f$u\f$是未知标量, \f$a\f$是依赖于坐标的向量,
 * \f$ a=\omega\otimes r\f$, \f$r\f$为位置向量.
 *   - 二维情形: \f$ a=(-y,x)^T \f$;
 *   - 三维情形: \f$\omega\f$从输入文件读入.
 *
 * \f$ u\f$定义在网格中心, 通量\f$ a*u\f$定义在网格面心.
 * 时间离散采用显式格式, 空间离散采用Goundov格式, 格式的精度由输入文件指定.
 *
 * 必须实现如下纯虚函数:
 *  - initializeComponent() : 初始化积分构件.
 *  - computeOnPatch() : 为数值构件完成数值计算.
 *  - getPatchDt()   : 为步长构件计算时间步长.
 *  - initializePatchData() : 为初值构件完成初始化.
 *
 * 细化数据通信时, 如果需要为数据片填充物理边界条件, 则必须重载如下函数:
 *  - setPhysicalBoundaryConditions() : 填充物理边界条件.
 *
 * 如果调用自定义细化插值算子或预(后)处理标准细化插值算子,
 * 则必须重载如下虚函数:
 *  - preprocessRefineOperator()  : 预处理细化插值算子. 
 *  - postprocessRefineOperator() : 后处理细化插值算子.
 *
 * 如果调用自定义粗化插值算子或预(后)处理标准粗化插值算子,
 * 则必须重载如下虚函数:
 *  - preprocessCoarsenOperator()  : 预处理粗化插值算子. 
 *  - postprocessCoarsenOperator() : 后处理粗化插值算子.
 *
 * @see algs::IntegratorComponent
 */

class RotAdv :
   public algs::StandardComponentPatchStrategy<NDIM>,
   public tbox::Serializable
{
public:

enum ROTATION_MODEL {
                   CIRCLE_ROTATION_DELTA,
                   CIRCLE_ROTATION_PARAB,
                   SQUARE_ROTATION_DELTA
                 };

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
   RotAdv(const string& object_name,
          tbox::Pointer<tbox::Database> input_db,
          tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom);

   /*! @brief 虚构函数.
    */
   virtual ~RotAdv();

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
    *
    */
   virtual void initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name);

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

   // 计算通量.
   void computeFluxOnPatch(hier::Patch<NDIM>& patch,
                           const double time,
                           const double dt);

   // 守恒型差分, 由通量校正数值解.
   void conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const double dt);

   /*!
    * @brief 标记待细化的网格单元.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param tag_time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param tag_index 输入参数, 整型, 表示存储标记值的数据片.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    */
   virtual void tagCellsOnPatch(hier::Patch<NDIM>& patch,
                                const double  tag_time,
                                const int     tag_index,
                                const bool    initial_time);

#if (NDIM==3)
   /*
    * 两种三维校正通量的计算方法. 由函数updateFluxesOnPatch()调用.
    */
   void compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                            const double dt);
   void compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                            const double dt);
#endif


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
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > d_javis_writer;

   /**
    * 存储数值解 - [u]
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_uval;
   /**
    * 存储通量  - [F=a*u]
    */
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_flux;

   // 上下文: 当前值, 新值, 可视化.
   tbox::Pointer< hier::VariableContext > d_current;
   tbox::Pointer< hier::VariableContext > d_new;
   tbox::Pointer< hier::VariableContext > d_scratch;
   tbox::Pointer< hier::VariableContext > d_plot_context;

   // 数据片索引号. 
   int d_uval_current_id, d_uval_new_id, d_uval_scratch_id, d_flux_new_id;

   /*
    *  用户数值子程序中使用的参数:
    *    d_godunov_order ....... Godunov格式的阶数(1, 2, 或4)
    *    d_corner_transport .... 三维有限差分中的通量校正方式.
    *    d_nghosts ............. 中心量和面心量变量的影像区宽度(为单元为单位).
    *    d_fluxghosts .......... 通量(面心量)变量的影像区宽度(为单元为单位).
    */
   int d_godunov_order;
   string d_corner_transport;
   hier::IntVector<NDIM> d_nghosts;
   hier::IntVector<NDIM> d_fluxghosts;

#if (NDIM==3)
   /*
    * 刚体旋转向量w, 由输入文件指定. 该参数仅适应三维.
    */
   double d_omega[NDIM];
#endif

   /*
    * 物理模型.
    *    CIRCLE_ROTATION_DELTA: 以d_center为中心, 以d_radius为半径的园,
    *                           园内值为1, 园外值为0.
    *    CIRCLE_ROTATION_PARAB: 以d_center为中心, 以d_radius为半径的园,
    *                           园内值为已知的抛物函数, 园外值为0.
    *    SQARE_ROTATION_DELTA : 以d_center为中心, 以d_radius为内接园半径的正方体,
    *                           园内值为已知的抛物函数, 园外值为0.
    */
   string d_model_problem;
   double d_center[NDIM];
   double d_radius;

   /*
    * 存储误差评估参数: 误差阈值.
    */
   tbox::Array<double> d_grad_tol;  // 梯度检测 :  |grad(u)|>d_grad_tol;

};

#endif
