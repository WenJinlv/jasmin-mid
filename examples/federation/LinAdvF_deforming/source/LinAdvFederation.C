//
// 文件名: LinAdvFederation.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 $
// 描述  : 线性对流问题的网格片时间积分算法类的实现.
//

#include "LinAdvFederation.h"

#include "CellData.h"
#include "NodeData.h"
#include "SideData.h"
#include "DeformingPatchGeometry.h"
#include "tbox/RestartManager.h"
using namespace JASMIN;

#include <assert.h>
using namespace std;


/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
LinAdvFederation::LinAdvFederation(const string& object_name,
               tbox::Pointer< hier::FederationGridGeometry<NDIM> > grid_geometry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!grid_geometry.isNull());
#endif

   d_object_name = object_name;
   d_grid_geometry = grid_geometry;

   d_federal_c   = -1;
   d_federal_f = -1;

   d_uval_current_id = -1;
   d_uval_scratch_id = -1;
   d_coords_id = -1;

   // 注册变量和数据片.
   registerModelVariables();

}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
LinAdvFederation::~LinAdvFederation()
{
tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};


/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void LinAdvFederation::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 上下文, 变量: 注册到数据库.
   d_current = variable_db->getContext("CURRENT");
   d_scratch = variable_db->getContext("SCRATCH");
   d_uval = variable_db->getVariable("uval");
   if(d_uval.isNull()) d_uval = new pdat::CellVariable<NDIM,double>("uval",1);
   d_uval_current_id =
        variable_db->registerVariableAndContext(d_uval, d_current);
   d_uval_scratch_id =
        variable_db->registerVariableAndContext(d_uval, d_scratch, 1);

   d_coords = variable_db->getVariable("coords");
   if(d_coords.isNull()) d_coords = new pdat::NodeVariable<NDIM,double>("coords",NDIM);
   d_coords_id =
        variable_db->registerVariableAndContext(d_coords, d_current);

   // 邦员的索引号.
   d_federal_c   = d_grid_geometry->getFederalNumber("LinAdv_C");
   d_federal_f   = d_grid_geometry->getFederalNumber("LinAdv_F");

}


/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void LinAdvFederation::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string &intc_name = intc->getName();

   if ( intc_name=="FEDERATION_COUPLE") {  // 联邦积分构件.
         intc->registerMeshMappingMethod(d_federal_c,
                                         d_federal_f,
                                         "UniCoarsen",
                                         d_coords_id);
         intc->registerCouplePatchDataF2C(d_uval_current_id,
                                       d_uval_current_id,
                                       "LINEAR_COUPLE",
                                       d_uval_scratch_id);   
         intc->registerCouplePatchDataC2F(d_uval_current_id,
                                       d_uval_current_id,
                                       "LINEAR_COUPLE",
                                       d_uval_scratch_id);   
    } else {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
    }

}

