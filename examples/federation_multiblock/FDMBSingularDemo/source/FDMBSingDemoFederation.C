//
// 文件名: FDMBSingDemoFederation.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 $
// 描述  : 验证联邦多块变形结构网格功能的网格片时间积分算法类.
//

#include "FDMBSingDemoFederation.h"

#include "CellData.h"
#include "SideData.h"
#include "tbox/RestartManager.h"
using namespace JASMIN;

#include <assert.h>
using namespace std;


/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
FDMBSingDemoFederation::FDMBSingDemoFederation(const string& object_name,
      tbox::Pointer< hier::FederationMultiblockGridGeometry<NDIM> > grid_geometry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!grid_geometry.isNull());
#endif

   d_object_name = object_name;
   d_grid_geometry = grid_geometry;

   d_federal_left   = -1;
   d_federal_middle = -1;
   d_federal_right  = -1;

   d_coords_current_id = -1;
   d_coords_scratch_id = -1;
   d_center_var_current_id = -1;
   d_center_var_scratch_id = -1;

   // 注册变量和数据片.
   registerModelVariables();

}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
FDMBSingDemoFederation::~FDMBSingDemoFederation()
{
tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};


/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void FDMBSingDemoFederation::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 上下文, 变量: 注册到数据库.
   d_current = variable_db->getContext("CURRENT");
   d_scratch = variable_db->getContext("SCRATCH");

   d_coords = variable_db->getVariable("coordinates");
   if(d_coords.isNull()) 
      d_coords = new pdat::NodeVariable<NDIM,double>("coordinates",NDIM);
   d_coords_current_id =
      variable_db->registerVariableAndContext(d_coords, d_current);
   d_coords_scratch_id =
      variable_db->registerVariableAndContext(d_coords, d_scratch, 2);
    
   d_center_var = variable_db->getVariable("center_var");
   if(d_center_var.isNull()) 
      d_center_var = new pdat::CellVariable<NDIM,double>("center_var",1);
   d_center_var_current_id =
      variable_db->registerVariableAndContext(d_center_var, d_current);
   d_center_var_scratch_id =
      variable_db->registerVariableAndContext(d_center_var, d_scratch, 1);

   // 邦员的索引号.
   d_federal_left   = d_grid_geometry->getFederalNumber("FDMBSingDemo_LEFT");
   d_federal_middle = d_grid_geometry->getFederalNumber("FDMBSingDemo_MIDDLE");
   d_federal_right  = d_grid_geometry->getFederalNumber("FDMBSingDemo_RIGHT");

}


/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void FDMBSingDemoFederation::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string &intc_name = intc->getName();

   if ( intc_name=="FEDERATION") {  // 联邦积分构件.
         intc->registerFederals(d_federal_left,  d_federal_middle);
         intc->registerFederals(d_federal_left,  d_federal_right);
         intc->registerFederals(d_federal_middle,d_federal_left);
         intc->registerFederals(d_federal_middle,d_federal_right);
         intc->registerFederals(d_federal_right, d_federal_left);
         intc->registerFederals(d_federal_right, d_federal_middle);
         intc->registerRefinePatchData(d_coords_scratch_id,
                                       d_coords_current_id);
         intc->registerRefinePatchData(d_center_var_scratch_id,
                                       d_center_var_current_id);
    } else {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
    }

}

