//
// 文件名:	LinAdvTagGeometry.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 250 $
// 修改  :	$Date: 2007-05-24 10:56:58$
// 描述  :	网格几何标记类
//

#ifndef included_LinAdvTagGeometry
#define included_LinAdvTagGeometry
 
#include "TagGeometryLevelStrategy.h"
#include "NumericalIntegratorComponent.h"
using namespace JASMIN;

#include "LinAdv.h"

/**
 * @brief 该类从网格几何标记策略类 
 * appu::TagGeometryLevelStrategy 派生,
 * 标记需要计算的索引范围，重建网格几何.
 */
class LinAdvTagGeometry 
: public appu::TagGeometryLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    * @param object_name          输入参数, 字符串, 表示对象名称.
    * @param patch_strategy       输入参数, 指针, 线性对流方程网格片积分算法.
    */     
   LinAdvTagGeometry(const string& object_name,
                        LinAdv* patch_strategy);

   /**
    * @brief 析构函数.
    */
   virtual ~LinAdvTagGeometry();

   ///@name 重载基类appu::TagGeometryLevelStrategy<NDIM>的函数
   //@{
   
   /**
    * @brief 初始化: 重建网格几何需要的积分构件.
    * @param manager 输入参数, 指针, 指向积分构件管理器.
    */
   void initializeLevelIntegrator(
           tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 标记需要保留的网格单元. 
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param tag_index  输入参数, 整型, 存储标记值的网格数据片的索引号.
    */
   void tagCellsForGeometry(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                const int tag_index);

   //@}

private:
   /*!@brief 对象名称. */
   string d_object_name; 

   /*!@brief 线性对流的网格片积分算法. */
   LinAdv* d_patch_strategy;

   /*!@brief 数值构件: 标记待细化的网格单元. */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_tag_intc;

};

#endif
