//
// 文件名: MBLinAdv.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 
// 描述  : 线性对流问题的网格片时间积分算法类.
//

#ifndef included_MBLinAdv
#define included_MBLinAdv

#include "StandardComponentPatchStrategy.h"
#include "MultiblockUniRectangularGridGeometry.h"
#include "CellVariable.h"
#include "SideVariable.h"
#include "JaVisDataWriter.h"
using namespace JASMIN;

/**
 * @brief 该类从标准构件网格片策略类 algs::StandardComponentPatchStrategy 派生, 
 * 实现求解线性对流方程的数值计算子程序.
 *
 * 该类需要从输入文件读取如下参数. 
 *     - \b constant_x_velocity \n
 *       双精度浮点型, X方向的对流速度.
 */
class MBLinAdv
: public algs::StandardComponentPatchStrategy<NDIM>,
  public tbox::Serializable
{
public:

   /*! @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示对象名称.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    * @param grid_geom   输入参数, 指针, 指向均匀矩形网格几何.
    *
    * @note
    * 该函数主要完成以下操作:
    *  -# 初始化内部数据成员;
    *  -# 定义变量和数据片.
    */
   MBLinAdv(const string& object_name,
          tbox::Pointer<tbox::Database> input_db,
          tbox::Pointer<geom::MultiblockUniRectangularGridGeometry<NDIM> > grid_geom);

   /*!
    * @brief 析构函数.
    */
   virtual ~MBLinAdv();

   /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
   // @{

   /*!
    * @brief 初始化指定的积分构件.
    *
    * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
    *
    * @param intc 输入参数, 指针, 指向待初始化的积分构件对象.
    */
    void initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const;

   /**
    * @brief 初始化数据片（支持初值构件）.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param intc_name 输入参数, 字符串, 当前调用该函数的初值构件之名称.
    */
   void initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name);

   /*!
    * @brief 完成单个网格片上的数值计算（支持数值构件）.
    *
    * 该函数基于显式迎风格式，实现通量和守恒量的计算。
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示当前时刻.
    * @param dt    输入参数, 双精度浮点型, 表示时间步长.
    * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
    * @param intc_name 输入参数, 字符串, 当前调用该函数的数值构件之名称.
    */
   void computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name);

   /*! 
    * @brief 在单个网格片上计算稳定性时间步长（支持时间步长构件）.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt   输入参数, 整型, 表示前次积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示前次积分的时间步长.
    * @param intc_name 输入参数, 字符串, 当前调用该函数的时间步长构件之名称.
    *
    * @return 双精度浮点型, 时间步长.
    */
   double getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name);

   /**
    * @brief 根据边界条件, 填充物理边界影像区数据.
    *
    * @param patch     输入参数, 网格片类, 表示网格片.
    * @param fill_time 输入参数, 双精度浮点型, 表示当前时刻.
    * @param ghost_width_to_fill 输入参数, 整型向量, 表示物理边界影像区被填充的宽度.
    * @param intc_name 输入参数, 字符串, 当前调用该函数的积分构件之名称.
    */
   void setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);

   //@}
   
   ///@name 重载基类tbox::Serializable的函数
   //@{
   /*!
    * @brief 将数据成员输出到重启动数据库.
    * @param db 输入参数, 指针, 指向重启动数据库.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   //@}
   
   ///@name 自定义函数
   //@{

   /*!
    * @brief 注册绘图量.
    * @param javis_writer 输入参数, 指针, 表示 JaVis 数据输出器.
    */
   void registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer);

   //@}

private:

   /*!@brief 从输入数据库读入数据.  */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   /*!@brief 从重启动数据库读入数据.  */
   void getFromRestart();
  
   /*!@brief 注册变量和数据片.  */
   void registerModelVariables();

   /*!@brief 对象名.  */
   string d_object_name;

   /*!@brief 网格几何. */
   tbox::Pointer<geom::MultiblockUniRectangularGridGeometry<NDIM> > d_grid_geometry;

   /*!@brief 守恒量 u */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_uval;

   /*!@brief 通量 f */
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_flux;

   /*!@brief 当前值上下文 */
   tbox::Pointer< hier::VariableContext > d_current;

   /*!@brief 新值上下文 */
   tbox::Pointer< hier::VariableContext > d_new;

   /*!@brief 演算值上下文 */
   tbox::Pointer< hier::VariableContext > d_scratch;

   /*!@brief <uval, current> 数据片的索引号 */
   int d_uval_current_id;

   /*!@brief <uval, new> 数据片的索引号 */
   int d_uval_new_id;

   /*!@brief <flux, new> 数据片的索引号 */
   int d_flux_new_id;

   /*!@brief <uval, scratch> 数据片的索引号 */
   int d_uval_scratch_id;

   /*!@brief x向速度常量 */
   double d_x_velocity;

};

#endif
