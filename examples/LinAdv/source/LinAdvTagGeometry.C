 //
// 文件名:  LinAdvTagGeometry.C
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  重构网格几何.
//


#include "LinAdvTagGeometry.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
LinAdvTagGeometry::LinAdvTagGeometry(
        const string& object_name,
        LinAdv* patch_strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(patch_strategy!=NULL);
#endif

   d_object_name = object_name;
   d_patch_strategy = patch_strategy;
}

/*************************************************************************
 *
 * 析构函数.
 *
 ************************************************************************/
LinAdvTagGeometry::~LinAdvTagGeometry() {

}

/*************************************************************************
 * 初始化: 创建所有计算需要的积分构件.
 *
 * 该函数创建了1个初值构件
 *
 *************************************************************************/
void LinAdvTagGeometry::initializeLevelIntegrator(
            tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

    //数值构件: 标记待细化的网格单元.
    d_tag_intc  = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "TAG_CELLS_GEOM",
                                                        d_patch_strategy,
                                                        manager);

}

/*************************************************************************
 *
 * 标记待细化的网格单元。
 *
 ************************************************************************/
void LinAdvTagGeometry::tagCellsForGeometry(
               const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
               const int tag_index)
{

   // 设置存储标记值的数据片索引号.
   d_patch_strategy->setTagIndex(tag_index);

   // 标记细化单元.
   d_tag_intc->computing(level,
                         0.0,
                         0.0);

   d_patch_strategy->clearTagIndex();

}

