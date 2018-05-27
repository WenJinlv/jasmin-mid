//
// 文件名:	MBEulerMM.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 170 $
// 修改  :	$Date: 2007-06-27 08:44:22 $
// 描述  :	Euler方程组：单网格片数值计算子程序.
//
 
#ifndef included_MBEulerMM
#define included_MBEulerMM

#ifndef included_JASMIN_config
#include "JASMIN_config.h"
#endif

#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Serializable.h"

#include "CellVariable.h"
#include "NodeVariable.h"
#include "CellData.h"
#include "NodeData.h"
#include "SideVariable.h"
#include "SideData.h"
#include "VariableContext.h"
#include "Patch.h"

#include "JaVisDerivedDataStrategy.h"
#include "JaVisDataWriter.h"
#include "DomainDataWriter.h"

#if (NDIM == 2)
#include "DeformingGridInputUtilities2.h"
#endif

#include "MultiblockDeformingGridGeometry.h"
#include "StandardComponentPatchStrategy.h"
#include "IntegratorComponent.h"

using namespace JASMIN;


/**
 * @brief 类 MBEulerMM 在单个网格片, 提供数值计算子程序, 
 * 基于移动网格方法求解Euler方程组.
 * 
 * 实现网格片时间积分算法基类 algs::StandardComponentPatchStrategy.
 * 在变形的单个网格片上, 求解Euler方程组.
 *
 * 实现如下纯虚函数:
 *  - initializeComponent() : 支撑所有积分构件的初始化.
 *  - computeOnpatch() :      支撑数值构件, 完成通信和计算.
 *  - getPatchDt()   :        支撑步长构件, 计算时间步长.
 *  - initializePatchData() : 为初值构件完成初始化.
 *
 * 填充物理边界条件, 则重载如下函数:
 *  - setPhysicalBoundaryConditions() : 填充物理边界条件.
 *
 * 由于移动网格方法仅考虑单层网格, 因此, 不需要插值算子.
 * 有关移动网格方法, 请参考主控程序的说明.
 *
 * @see algs::IntegratorComponent
 */

class MBEulerMM : 
   public tbox::Serializable,
   public algs::StandardComponentPatchStrategy<NDIM>
{
public:

   /*!
    * @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示该对象的名字.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    * @param grid_tool   输入参数, 指针, 指向变形网格辅助输入工具.
    * @param grid_geom   输入参数, 指针, 指向多块均匀矩形网格几何.
    * @note
    * 该函数执行以下操作:
    *  - 创建变量;
    *  - 初始化私有数据成员;
    *  - 注册对象到重启动管理器;
    *  - 如果需断点续算, 则从重启动文件中读入数据; 
    *  - 从输入数据库中读入数据, 该操作将覆盖从重启动文件中读入的数据.
    * 
    */

   MBEulerMM(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
         tbox::Pointer< geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom);  

    /**
     * @brief 析构函数.
     */
   ~MBEulerMM();


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

   virtual void reduceOnPatch( 
                    double *error,
                    int     len,
                    hier::Patch<NDIM> &patch,
                    const double time,
                    const double dt,
                    const string& intc_name );
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
   void registerDataWriter(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer,
      tbox::Pointer<appu::DomainDataWriter<NDIM> > domain_writer);

   /*!
    * @brief 输出当前对象中所有数据成员到指定输出流.
    * @param os 输入输出参数, 输出流类型, 表示指定的输出流.
    */
   void printClassData(ostream& os) const;

private:
 
   /*
    * 从输入数据库和重启动数据库读入数据. 
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);
   void getFromRestart();

   // 创建所有变量和上下文, 注册变量上下文配对到变量数据库.
   void registerModelVariables();

   /// @brief 以下函数在时间积分过程中被依次调用.
   //@{

   // 为监控函数赋初值.
   void computeWCNEW(hier::Patch<NDIM>& patch,
                     const double  time,
                     const double  dt,
                     const bool    initial_time);

   // 迭代使监控函数收敛.
   void wcnewIteration(hier::Patch<NDIM>& patch,
                       const double  time,
                       const double  dt,
                       const bool    initial_time);

   // 移动网格.
   void moveMesh(hier::Patch<NDIM>& patch,
                       const double  time,
                       const double  dt,
                       const bool    initial_time);

   // 网格移动后, 修正数值解.
   void moveSolution(hier::Patch<NDIM>& patch,
                       const double  time,
                       const double  dt,
                       const bool    initial_time);

   // 求方向导数.
   void computeGradient(hier::Patch<NDIM>& patch,
                        const double  time,
                        const double  dt,
                        const bool    initial_time);


   // 求Euler方程.
   void solveEulerEquation(hier::Patch<NDIM>& patch,
                           const double  time,
                           const double  dt,
                           const bool    initial_time);

   //@}

   // 以下为私有数据成员.

   /*
    * 对象名.
    */
   string d_object_name;     

   // 变形网格辅助输入工具.
   tbox::Pointer< appu::DeformingGridInputUtilities2 > d_grid_tool;

   /*
    * 网格几何、可视化数据输出器和绘图上下文 .
    */
   tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > d_grid_geometry;
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > d_javis_writer;

   // Euler方程的参数: 理想气体比热比.
   double d_gamma;

   // 移动网格的控制参数
   double d_alpha, d_beta;

   // 计算区域的参数.
   double d_xlo[NDIM],d_xhi[NDIM];
   hier::IntVector<NDIM> d_ncells;

   // 存储数值解: 坐标,密度,动量,能量密度.
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_coords;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_solution;

   // 数值解辅助量: 通量, 移动网格的控制函数.
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_flux;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_wc;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_dux;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_duy;

   // 上下文: 当前值, 新值, 演算值.
   tbox::Pointer< hier::VariableContext > d_current;
   tbox::Pointer< hier::VariableContext > d_new;
   tbox::Pointer< hier::VariableContext > d_scratch;
   tbox::Pointer< hier::VariableContext > d_move;

   // 数据片索引号.
   int d_coords_current_id,   // < d_coords, d_current >: 存储结点坐标的当前值.
       d_coords_new_id,       // < d_coords, d_new >    : 存储结点坐标的最新值.
       d_coords_scratch_id,   // < d_coords, d_scratch >: 存储结点坐标的演算值.
       d_coords_move_id;      // < d_coords, d_move >   : 存储结点坐标的移动值.
   int d_solution_current_id,  // < d_solution, d_current >:存储数值解的当前值.
       d_solution_new_id,      // < d_solution, d_new > : 存储数值解的最新值.
       d_solution_scratch_id;  // < d_solution, d_scratch > :存储数值解的演算值.

   int d_flux_new_id;         // < d_flux, d_new > : 存储通量.
   int d_wc_new_id,           // < d_wc, d_new >   : 存储控制函数值的最新值.
       d_wc_scratch_id;       // < d_wc, d_scratch > : 存储控制函数值的演算值.
   int d_dux_scratch_id,      // < d_dux, d_scratch >: 存储X方向导数值. 
       d_duy_scratch_id;      // < d_duy, d_scratch >: 存储Y方向导数值.

   // 数据片的两种影像区宽度.
   hier::IntVector<NDIM> d_zeroghosts, d_oneghosts;
                  
};

#endif
