//
// 文件名: LinAdvFederation.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 
// 描述  : 线性对流问题的网格片时间积分算法类.
//

#ifndef included_LinAdvFederation
#define included_LinAdvFederation

#include "StandardComponentPatchStrategy.h"
#include "FederationGridGeometry.h"
#include "CellVariable.h"
#include "NodeVariable.h"
using namespace JASMIN;

/**
 * @brief 该类从标准构件网格片策略类 algs::StandardComponentPatchStrategy 派生, 
 * 为联邦计算实现网格片策略类.
 */
class LinAdvFederation
: public algs::StandardComponentPatchStrategy<NDIM>
{
public:

   /*! @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示对象名称.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    * @param grid_geom   输入参数, 指针, 指向联邦均匀矩形网格几何.
    *
    * @note
    * 该函数主要完成以下操作:
    *  -# 初始化内部数据成员;
    *  -# 定义变量和数据片.
    */
   LinAdvFederation(const string& object_name,
          tbox::Pointer< hier::FederationGridGeometry<NDIM> > grid_geometry);

   /*!
    * @brief 析构函数.
    */
   virtual ~LinAdvFederation();

   /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
   // @{

   /*!
    * @brief 初始化联邦积分构件.
    *
    * @param intc 输入参数, 指针, 指向待初始化的联邦积分构件对象.
    */
    void initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const;

private:

   /*!@brief 注册变量和数据片.  */
   void registerModelVariables();

   /*!@brief 对象名.  */
   string d_object_name;

   /*!@brief 联邦网格几何. */
   tbox::Pointer< hier::FederationGridGeometry<NDIM> > d_grid_geometry;

   /*!@brief 守恒量 u */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_uval;
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_coords;

   /*!@brief 当前值上下文 */
   tbox::Pointer< hier::VariableContext > d_current;

   /*!@brief 演算值上下文 */
   tbox::Pointer< hier::VariableContext > d_scratch;

   /*!@brief <uval, current> 数据片的索引号 */
   int d_uval_current_id;
   int d_coords_id;

   /*!@brief <uval, scratch> 数据片的索引号 */
   int d_uval_scratch_id;

   /*!@brief 邦员的索引号. */
   int d_federal_c, d_federal_f;

};

#endif
