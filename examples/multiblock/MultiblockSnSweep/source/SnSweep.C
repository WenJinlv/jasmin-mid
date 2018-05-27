//
// 文件名: SnSweep.C
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 $  
// 描述  : Sn输运计算的网格片策略类实现.
//

#include "SnSweep.h"
#include "SnSweepFort.h"

#include "VariableDatabase.h"
#include "GroupNodeVariable.h"
#include "GroupNodeData.h"
#include "NodeVariable.h"
#include "NodeData.h"
#include "SideData.h"
#include "BlockUniRectangularPatchGeometry.h"

#include "tbox/TimerManager.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// 定义宏
#if (NDIM==2)
#define BOX_INDEX(BOX)                                \
        BOX.lower()(0),BOX.upper()(0),                \
        BOX.lower()(1),BOX.upper()(1)
#define ABC_DRESS                                     \
        am.getPointer(), bm.getPointer()
#elif (NDIM==3)
#define BOX_INDEX(BOX)                                \
        BOX.lower()(0),BOX.upper()(0),                \
        BOX.lower()(1),BOX.upper()(1),                \
        BOX.lower()(2),BOX.upper()(2)   
#define ABC_DRESS                                     \
        am.getPointer(), bm.getPointer(), cm.getPointer()
#endif


/****************************************************************
* 构造函数.
*****************************************************************/
SnSweep::SnSweep(const string& object_name,
               tbox::Pointer<tbox::Database> input_db,
               tbox::Pointer<geom::MultiblockUniRectangularGridGeometry<NDIM> > grid_geom)
: d_object_name(object_name),
  d_grid_geometry(grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
#endif

   // 从输入文件读入参数
   getFromInput(input_db);

   // 配置离散方向
   d_angles.resizeArray(d_ms*NDIM);
   d_sn_weight.resizeArray(d_ms);
   snquad_(d_angles.getPointer(),
           d_sn_weight.getPointer(),
           d_sn_cos.getPointer(),
           d_ms,d_sn_n);
   
   // 创建所有变量及数据片索引号
   registerModelVariables();

   // 
   tbox::Pointer< geom::UniRectangularGridGeometry<NDIM> > bgeom = 
                  grid_geom->getBlockGeometry(0);
   const double *xup = bgeom->getXUpper();
   const double *xlo = bgeom->getXLower();
   for ( int i=0; i<NDIM; i++ ) {
      d_xleng[i] = xup[i] - xlo[i]; 
   }

}

/****************************************************************
* 构造函数.
*****************************************************************/
SnSweep::~SnSweep()
{
};

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/
void SnSweep::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   // 1. fi    : 角通量.
   // 2. sgm_t : 输运截面
   // 3. sgm_s : 散射截面
   // 4. qs    : 散射源
   // 6. qnw   : 数值外源项
   // 7. wf    : 流通量
   tbox::Pointer<hier::Variable<NDIM> > var_fi
              = new pdat::GroupNodeVariable<NDIM,double>("fi",d_mg,d_ms); // GroupNode类型

   tbox::Pointer<hier::Variable<NDIM> > var_sgm_t
              = new pdat::GroupNodeVariable<NDIM,double>("sgm_t",d_mg,1); // GroupNode类型

   tbox::Pointer<hier::Variable<NDIM> > var_sgm_s
              = new pdat::GroupNodeVariable<NDIM,double>("sgm_s",d_mg*(d_mg+1)/2,1); // GroupNode类型

   tbox::Pointer<hier::Variable<NDIM> > var_qs
              = new pdat::GroupNodeVariable<NDIM,double>("qs",d_mg, 1); // GroupNode类型

   tbox::Pointer<hier::Variable<NDIM> > var_qnw
              = new pdat::GroupNodeVariable<NDIM,double>("qnw",d_mg, d_ms); // GroupNode类型

   tbox::Pointer<hier::Variable<NDIM> > var_wf
              = new pdat::NodeVariable<NDIM,double>("wf",d_mg); // Node类型

   // 定义上下文.
   tbox::Pointer<hier::VariableContext> context_cur
                        = variable_db -> getContext("CURRENT");
   tbox::Pointer<hier::VariableContext> context_new
                        = variable_db -> getContext("NEW");
   tbox::Pointer<hier::VariableContext> context_scr
                        = variable_db -> getContext("SCRATCH");


   // 申请数据片索引号.
   d_fi_new_id = variable_db->registerVariableAndContext(var_fi,context_new,0);
   d_fi_cur_id = variable_db->registerVariableAndContext(var_fi,context_cur,0);
   d_sgm_t_id  = variable_db->registerVariableAndContext(var_sgm_t,context_cur,0);
   d_sgm_s_id  = variable_db->registerVariableAndContext(var_sgm_s,context_cur,0);
   d_qs_scr_id = variable_db->registerVariableAndContext(var_qs,context_scr,0);
   d_qnw_scr_id = variable_db->registerVariableAndContext(var_qnw,context_scr,0);
   d_wf_cur_id = variable_db->registerVariableAndContext(var_wf,context_cur,0);
   d_wf_new_id = variable_db->registerVariableAndContext(var_wf,context_new,0);
   d_wf_scr_id = variable_db->registerVariableAndContext(var_wf,context_scr,0);

}

/********************************************************************
*  初始化积分构件.
*********************************************************************/
void SnSweep::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件. 
        intc->registerInitPatchData(d_fi_cur_id);
        intc->registerInitPatchData(d_wf_cur_id);
        intc->registerInitPatchData(d_sgm_t_id);
        intc->registerInitPatchData(d_sgm_s_id);
                                              
   }else if(intc_name=="SOLVE") { // 扫描构件: 隐式迎风格式求解角通量.
        intc->registerPatchData(d_fi_new_id);

   }else if(intc_name=="SOURCEW") { // 数值构件: 计算外源项，同时将复制current-->new
        intc->registerRefinePatchData(d_fi_new_id,
                                    d_fi_cur_id);
        intc->registerRefinePatchData(d_wf_new_id,
                                    d_wf_cur_id);

   }else if(intc_name=="NEW_2_CUR") { // new-->current的复制构件
        intc->registerCopyPatchData(d_fi_cur_id,
                                    d_fi_new_id);
        intc->registerCopyPatchData(d_wf_cur_id,
                                    d_wf_new_id);

   }else if(intc_name=="NEW"){// 为新时刻未知量数据片调度内存.
        intc->registerPatchData(d_fi_new_id);
        intc->registerPatchData(d_wf_new_id);

   }else if(intc_name=="SCRATCH"){// 为新时刻未知量数据片调度内存.
        intc->registerPatchData(d_qnw_scr_id);
        intc->registerPatchData(d_qs_scr_id);
        intc->registerPatchData(d_wf_scr_id);

   }else if( intc_name!="SCATTER" 
#ifdef TESTING
          && intc_name!="DIFF_RESULT"
#endif
          && intc_name!="VFLUX"
          && intc_name!="ITER_ERROR" ) {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
   }

}


/********************************************************************
*  初始化数据片.
*********************************************************************/
void SnSweep::initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT");
#endif

   if (initial_time) {

      tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_cur = patch.getPatchData(d_fi_cur_id);
      tbox::Pointer< pdat::NodeData<NDIM,double> >      wf_cur = patch.getPatchData(d_wf_cur_id);
      tbox::Pointer< pdat::GroupNodeData<NDIM,double> > sgm_t = patch.getPatchData(d_sgm_t_id);
      tbox::Pointer< pdat::GroupNodeData<NDIM,double> > sgm_s = patch.getPatchData(d_sgm_s_id);

      const hier::Box<NDIM>& patch_box = patch.getBox();
      const tbox::Pointer<geom::BlockUniRectangularPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();

      tbox::Array<double> am(2*(1+patch_box.numberCells(0)));
      tbox::Array<double> bm(2*(1+patch_box.numberCells(1)));
      #if NDIM==3
      tbox::Array<double> cm(2*(1+patch_box.numberCells(2)));
      #endif

      int mg1 = d_mg*(d_mg+1)/2;
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(sgm_t->getGhostCellWidth()==hier::IntVector<NDIM>(0));
      assert(sgm_s->getGhostCellWidth()==hier::IntVector<NDIM>(0));
      assert(wf_cur->getGhostCellWidth()==hier::IntVector<NDIM>(0));
      assert(fi_cur->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
      init_( BOX_INDEX(patch_box),
             d_mg, mg1,
             d_ms,
             sgm_t->getPointer(),
             sgm_s->getPointer(), 
             wf_cur->getPointer(),    
             fi_cur->getPointer(),   
             dx, d_xleng,
             d_sgm_t0.getPointer(), d_sgm_s0.getPointer(),
             d_vn.getPointer(),
             d_angles.getPointer(),
             d_sn_weight.getPointer(),
             ABC_DRESS);

   } 
}

/********************************************************************
 *  为数值构件完成数值计算.
 ********************************************************************/
void SnSweep::computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name)
{
   if(intc_name=="SCATTER") {  // 计算散射源项.
       computeScatterItemOnPatch(patch); 

   } else if(intc_name=="SOURCEW") { // 计算外源项.
       computeSoureceItemOnPatch(patch,dt); 

   } else if(intc_name=="VFLUX") { // 更新流通量.
       updateVolumeFluxOnPatch(patch);

#ifdef TESTING
   } else if(intc_name=="DIFF_RESULT") { // 比对数值解和精确解.
       diffResultOnPatch(patch);
#endif

   } else {
       TBOX_ERROR("\n::computeOnComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
   }
}

/********************************************************************
 *  计算散射源
 ********************************************************************/
void SnSweep::computeScatterItemOnPatch(hier::Patch<NDIM>& patch)
{
   tbox::Pointer< pdat::NodeData<NDIM,double> >   wf = patch.getPatchData(d_wf_new_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > qs = patch.getPatchData(d_qs_scr_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > sgm_s = patch.getPatchData(d_sgm_s_id);

   const hier::Box<NDIM>& patch_box = patch.getBox();
   int mg1 = d_mg*(d_mg+1)/2;

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(wf->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(qs->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(sgm_s->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   scatter_(qs->getPointer(),
            sgm_s->getPointer(),
            wf->getPointer(),
            BOX_INDEX(patch_box),
            d_mg,
            mg1);

}

/********************************************************************
 *  计算外源项
 ********************************************************************/
void SnSweep::computeSoureceItemOnPatch(hier::Patch<NDIM>& patch, 
                                        const double dt)
{
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > qnw    = patch.getPatchData(d_qnw_scr_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_cur = patch.getPatchData(d_fi_cur_id);
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(qnw->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(fi_cur->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   #endif

   const hier::Box<NDIM>& patch_box = patch.getBox();

   // 设置外源项
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > qw 
        = new pdat::GroupNodeData<NDIM,double>(patch_box,d_mg,d_ms,0);
#ifndef TESTING
   qw->fillAll(0.e0);

#else
   const tbox::Pointer<geom::BlockUniRectangularPatchGeometry<NDIM> > pgeom 
                     = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();

   tbox::Array<double> am(2*(1+patch_box.numberCells(0)));
   tbox::Array<double> bm(2*(1+patch_box.numberCells(1)));
   #if NDIM==3
   tbox::Array<double> cm(2*(1+patch_box.numberCells(2)));
   #endif

   double yd = pow(d_yd1,-d_nsteps-1);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(qw->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   sourcew_(
            BOX_INDEX(patch_box),
            d_mg,
            d_ms,
            qw->getPointer(),
            dx,
            d_xleng,
            d_angles.getPointer(),
            d_vn.getPointer(),
            d_ysgm_s.getPointer(),
            d_ysgm_t.getPointer(),
            yd,
            ABC_DRESS);
#endif

   // 计算数值源项: qnw = fi/v_g/dt + qw
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(qw->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(qnw->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(fi_cur->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   sourceqnw_(
              BOX_INDEX(patch_box),
              d_mg,
              d_ms,
              dt,
              qnw->getPointer(),
              qw->getPointer(),
              fi_cur->getPointer(),
              d_vn.getPointer());

}


/********************************************************************
 * 为扫描构件计算单元间依赖关系.
 ********************************************************************/
void SnSweep::setCellsDataDependencyOnPatch(
                                JASMIN::pdat::SideData<NDIM, int>& cell_depend,
				hier::Patch<NDIM> &patch,
				const double time,
				const double dt,
				const string& intc_name )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="SOLVE");
#endif

   const hier::Box<NDIM>& patch_box = patch.getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cell_depend.getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   setdependency_(
           BOX_INDEX(patch_box),
           d_ms,
           cell_depend.getPointer(0),
           cell_depend.getPointer(1),
#if NDIM==3
           cell_depend.getPointer(2),
#endif
           d_angles.getPointer());

}


/********************************************************************
 *  实现扫描构件的网格片策略类成员函数,
 *  按隐式迎风格式计算角通量。
 ********************************************************************/
void SnSweep::sweepingOnPatch(hier::Patch<NDIM>& patch,
                               const int angle_id,
                               const int ncells,
                               const int cells[][NDIM],
                               const double time, 
                               const double dt, 
                               const JASMIN::pdat::SideData<NDIM, int>& cell_depend,
                               const string& intc_name )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="SOLVE");
#endif
#if 0
   tbox::Pointer<tbox::Timer> t_sweep_on_patch = tbox::TimerManager::getManager()->
                         getTimer("apps::SnSweep::sweepingOnPatch()");
   t_sweep_on_patch->start();
#endif

   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_new = patch.getPatchData(d_fi_new_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > qnw    = patch.getPatchData(d_qnw_scr_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > qs     = patch.getPatchData(d_qs_scr_id);
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > sgm_t  = patch.getPatchData(d_sgm_t_id);

   const hier::Box<NDIM>& patch_box = patch.getBox();
   const tbox::Pointer<geom::BlockUniRectangularPatchGeometry<NDIM> > pgeom = 
                            patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();


#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fi_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(qnw->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(qs->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(sgm_t->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   solve_( BOX_INDEX(patch_box),
           d_mg, d_ms,
           angle_id,
           d_angles.getPointer(),
           ncells,
           cells,
           fi_new->getPointer(),    
           qnw->getPointer(),
           qs->getPointer(),
           sgm_t->getPointer(),
           dt, dx,
           d_vn.getPointer() );

#if 0
   t_sweep_on_patch->stop();
#endif
}

/********************************************************************
 *  更新流通量.
 ********************************************************************/
void SnSweep::updateVolumeFluxOnPatch(hier::Patch<NDIM>& patch)
{
   const hier::Box<NDIM>& patch_box = patch.getBox();

   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_new = patch.getPatchData(d_fi_new_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> >      wf_new = patch.getPatchData(d_wf_new_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > wf_scr = patch.getPatchData(d_wf_scr_id);

   wf_scr->copy(*wf_new);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fi_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(wf_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   vflux_( BOX_INDEX(patch_box),
           d_mg,
           d_ms,
           wf_new->getPointer(),    
           fi_new->getPointer(),    
           d_sn_weight.getPointer() );
}

/********************************************************************
 * 计算流通量迭代误差.
 ********************************************************************/
void SnSweep::reduceOnPatch( 
                    double *vector,
                    int     len,
                    hier::Patch<NDIM> &patch,
                    const double time,
                    const double dt,
                    const string& intc_name )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="ITER_ERROR");
#endif
   const hier::Box<NDIM>& patch_box = patch.getBox();

   tbox::Pointer< pdat::NodeData<NDIM,double> > wf_new = patch.getPatchData(d_wf_new_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > wf_scr = patch.getPatchData(d_wf_scr_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(wf_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
   assert(wf_scr->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   vfluxerror_( 
           BOX_INDEX(patch_box),
           d_mg,
           wf_scr->getPointer(),    
           wf_new->getPointer(),
           vector);
}

/********************************************************************
 *  填充物理边界条件.
 ********************************************************************/
void SnSweep::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   if ( intc_name=="SOLVE" ) { // 填充边界条件

       tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_new = patch.getPatchData(d_fi_new_id);

       const tbox::Pointer<geom::BlockUniRectangularPatchGeometry<NDIM> > pgeom =patch.getPatchGeometry();
       const hier::Box<NDIM>& interior = patch.getBox();

       // 只处理棱型（二维）或面型（三维）物理边界.
       int btype = 1;
       hier::BoxArray<NDIM> fill_boxes;
       tbox::Array<int> loc_indexes;
       pgeom->getPhysicalBoundaryFillBoxes(
                                   fill_boxes,
                                   loc_indexes,
                                   btype, 
                                   interior,
				   hier::IntVector<NDIM>(1));
       for (int i = 0; i < fill_boxes.size(); i++) {
            fi_new->fillAll(0.0,fill_boxes(i));
       }
   }

}

/********************************************************************
 *  从输入文件读取数据成员. 
 ********************************************************************/
void SnSweep::getFromInput(tbox::Pointer<tbox::Database> db )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   d_matter_code= db->getString("matter");

   // 读入能群参数
   d_mg = db->getInteger("ngroup");
   d_vn.resizeArray(d_mg);
   d_vn = db->getDoubleArray("v");
   char matter_db_name[6+d_matter_code.size()];
   sprintf(matter_db_name,"Matter%s",d_matter_code.c_str());
   tbox::Pointer<tbox::Database> matter_db = db->getDatabase("Matter"+d_matter_code);
   d_sgm_t0.resizeArray(d_mg);
   d_sgm_t0 = matter_db->getDoubleArray("SGMt");
   int mg1 = d_mg*(d_mg+1)/2;
   d_sgm_s0.resizeArray(mg1);
   d_sgm_s0 = matter_db->getDoubleArray("SGMs");


   // 读入离散纵标系数
   d_sn_n = db->getInteger("s_n");
   int MAX_INT_BITS = 10;  //最大int型数的位数。
   char sn_db_name[1+MAX_INT_BITS+1];
   sprintf(sn_db_name,"S%d",d_sn_n);
   tbox::Pointer<tbox::Database> sn_db = db->getDatabase(sn_db_name);
   d_sn_cos.resizeArray(d_sn_n/2);
   d_sn_cos = sn_db->getDoubleArray("u");
   int ms1 = d_sn_n*(d_sn_n+2)/8;
   d_sn_weight.resizeArray(ms1);
   d_sn_weight = sn_db->getDoubleArray("w");
#if NDIM==2 
   d_ms = ms1 * 4;
#elif NDIM==3 
   d_ms = ms1 * 8;
#endif
}


#ifdef TESTING
/********************************************************************
 *  设置当前时间步数
 ********************************************************************/
 void SnSweep::setTimeStepNumber(const int nsteps)
{
   d_nsteps = nsteps; 
}

/********************************************************************
 *  设置精确解系数 yd0,yd1
 ********************************************************************/
 void SnSweep::setSolutionCoeff(const double dt)
{
   d_yd0 = 1000.;
   d_yd1= 1e0 + d_yd0 * dt;

   d_ysgm_t.resizeArray(d_mg);
   d_ysgm_s.resizeArray(d_mg);

   int lg = 0;
   for ( int mg=0; mg<d_mg; mg++) {
         d_ysgm_t[mg]=-d_yd0+d_sgm_t0[mg]*d_vn[mg];
         d_ysgm_s[mg]=0e0;
         for ( int ig=0; ig<=mg; ig++) {
            d_ysgm_s[mg]+=d_sgm_s0[lg]*d_vn[ig];
            lg++;
         }
#if NDIM==2
         d_ysgm_s[mg] *= 0.25;
#elif NDIM==3
         d_ysgm_s[mg] *= 0.125;
#endif
   }
}

/********************************************************************
 *  比对数值解和精确解
 ********************************************************************/
 void SnSweep::diffResultOnPatch(hier::Patch<NDIM> &patch)
{
   tbox::Pointer< pdat::GroupNodeData<NDIM,double> > fi_new = patch.getPatchData(d_fi_new_id);

   const hier::Box<NDIM>& patch_box = patch.getBox();
   const tbox::Pointer<geom::BlockUniRectangularPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();

   tbox::Array<double> am(2*(1+patch_box.numberCells(0)));
   tbox::Array<double> bm(2*(1+patch_box.numberCells(1)));
   #if NDIM==3
   tbox::Array<double> cm(2*(1+patch_box.numberCells(2)));
   #endif

   double yd = pow(d_yd1,-d_nsteps-1);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fi_new->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   diffsolution_(
            BOX_INDEX(patch_box),
            d_mg, d_ms,
            fi_new->getPointer(),
            dx, d_xleng,
            yd, d_vn.getPointer(), 
            d_angles.getPointer(),
            ABC_DRESS);
}

#endif

/*
*************************************************************************
*                                                                       *
* 注册绘图量到JaVis数据输出器.                                          *
*                                                                       *
*************************************************************************
*/
void SnSweep::registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(javis_writer.isNull()));
#endif

#ifdef HAVE_HDF5
   // 注册可视化量.
   javis_writer->registerDerivedPlotQuantity("Total Energy", "SCALAR", this, 1, 1.0, "NODE");
   javis_writer->registerPlotQuantity("FI", "COMPOSITE", d_fi_cur_id);
#endif
}

bool SnSweep::packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch<NDIM>& patch,
      const hier::Box<NDIM>& region,
      const string& variable_name,
      int depth_id)
{
   tbox::Pointer< pdat::NodeData<NDIM,double> >  wf = patch.getPatchData(d_wf_cur_id);

   const hier::Box<NDIM>& patch_box = patch.getBox();
   
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(wf->getGhostCellWidth()==hier::IntVector<NDIM>(0));
#endif
   computeenergy_( BOX_INDEX(patch_box),
           BOX_INDEX(region),
           d_mg, 
           wf->getPointer(),
           buffer );
   return true;
}



