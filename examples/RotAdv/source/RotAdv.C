//
// 文件名: RotAdv.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 (浜, 28  9 2007) $  
// 描述  : 刚体旋转问题的网格片时间积分算法类的实现.
//

#ifndef included_RotAdv_C
#define included_RotAdv_C

#include "RotAdv.h"

#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/Array.h"
#include "VariableDatabase.h"
#include "BoxArray.h"
#include "CellData.h"
#include "FaceData.h"
#include "VariableDatabase.h"

#include "UniRectangularPatchGeometry.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// Fortran数值子程序接口.
#include "RotAdvFort.h"

// Godunov格式精度, 变量的影像区宽度.
#define CELLG           (4)   // 中心量影像区宽度.
#define FACEG           (4)   // 面心量影像区宽度.
#define FLUXG           (1)   // 通量型变量影像区宽度.

// 逻辑值常量.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

using namespace std;

/****************************************************************
* 构造函数.
*****************************************************************/
RotAdv::RotAdv(const string& object_name,
               tbox::Pointer<tbox::Database> input_db,
               tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
#endif

   d_object_name = object_name;

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
   d_grid_geometry = grid_geom;

   // Godunov离散格式的缺省参数.
   d_godunov_order = 1;   // 格式的阶.
   d_corner_transport = "CORNER_TRANSPORT_1";
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(CELLG == FACEG);
#endif
   d_nghosts = hier::IntVector<NDIM>(CELLG);
   d_fluxghosts = hier::IntVector<NDIM>(FLUXG);

   // 读取输入参数和重启动数据, 并用其初始化对象.
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   // 传递常量给F77子程序.
   stufprobc_(CIRCLE_ROTATION_DELTA,CIRCLE_ROTATION_PARAB,
              SQUARE_ROTATION_DELTA,
              CELLG,FACEG,FLUXG);

   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();

}

/****************************************************************
* 析构函数.
*****************************************************************/
RotAdv::~RotAdv()
{
tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/
void RotAdv::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   d_uval    = new pdat::CellVariable<NDIM,double>("uval",1);  // 深度为1.
   d_flux    = new pdat::FaceVariable<NDIM,double>("flux",1);  // 深度为1.
  
   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");
   d_plot_context = d_current;

   const hier::IntVector<NDIM> zeroghosts(0),oneghosts(1);
   // 存储当前时刻的变量值: (uval,current), 影像区宽度为0.
   d_uval_current_id = variable_db->registerVariableAndContext(d_uval,
                                                           d_current,
                                                           zeroghosts);

   // 存储新时刻的变量值: (uval,new), 影像区宽度为格式宽度.
   d_uval_new_id = variable_db->registerVariableAndContext(d_uval,
                                                           d_new,
                                                           d_nghosts);
   // 存储变量的演算值: (uval,new), 影像区宽度为等于1.
   d_uval_scratch_id = variable_db->registerVariableAndContext(d_uval,
                                                           d_scratch,
                                                           oneghosts);
   // 存储新时刻的通量值: (flux,new), 影像区宽度等于1.
   d_flux_new_id = variable_db->registerVariableAndContext(d_flux,
                                                           d_new,
                                                           d_fluxghosts);

}

/*
*************************************************************************
*                                                                       *
* 注册JaVis数据输出器, 以采用JaVis工具对绘图文件进行后处理.           *
*                                                                       *
*************************************************************************
*/
void RotAdv::registerJaVisDataWriter(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(viz_writer.isNull()));
#endif
   d_javis_writer = viz_writer;

   // 注册可视化量.
   #ifdef HAVE_HDF5
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();
   d_javis_writer-> registerPlotQuantity("U",
                                         "SCALAR",
                       variable_db->mapVariableAndContextToIndex(
                                                d_uval,
                                                d_plot_context));
   #endif

}

/********************************************************************
*  初始化积分构件.
*********************************************************************/
void RotAdv::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件. 
        intc->registerInitPatchData(d_uval_current_id,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                    d_uval_scratch_id);
   }else if(intc_name=="FLUX") { // 计算通量的数值构件.
        intc->registerRefinePatchData(d_uval_new_id,
                                      d_uval_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE",
                                      d_uval_new_id,
                                      d_uval_current_id,
                                      d_uval_new_id);
   }else if(intc_name=="SYNC_FLUX") { // 沿粗细网格层交接面校正通量的同步构件. 
        intc->registerCoarsenPatchData(d_flux_new_id,
                                       "CONSERVATIVE_COARSEN");
   }else if(intc_name=="SYNC_OVERLAY") { // 在细网格层覆盖的内部区域, 同步.
        intc->registerCoarsenPatchData(d_uval_new_id,
                                       "CONSERVATIVE_COARSEN",
                                       d_uval_new_id);
   }else if(intc_name=="TAGC"){  // 标记构件.
        intc->registerRefinePatchData(d_uval_scratch_id,
                                      d_uval_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE");
   }else if(intc_name=="RESET_SOLUTION") { // 接收数值解.
        intc->registerCopyPatchData(d_uval_current_id,
                                    d_uval_new_id);
   }else if(intc_name=="SYNC_COPY") { // 同步前复制数据片.
        intc->registerCopyPatchData(d_uval_new_id,
                                    d_uval_current_id);
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 为新值数据片调度内存空间的构件.
        intc->registerPatchData(d_uval_new_id);
        intc->registerPatchData(d_flux_new_id);
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 为演算数据片调度内存空间的构件.
        intc->registerPatchData(d_uval_scratch_id);
   }else if(intc_name!= "TIME_STEP_SIZE" && 
            intc_name!="CONSERVATIVE_DIFFERENCE" ){
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
   }

}

/********************************************************************
*  计算时间步长.
*********************************************************************/
double RotAdv::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TIME_STEP_SIZE");
#endif
   
   NULL_USE(flag_last_dt);
   NULL_USE(last_dt);

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom
                       = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* xlo = pgeom->getXLower();

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
      patch.getPatchData(d_uval_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!uval_current.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = uval_current->getGhostCellWidth();

      double stabdt;
   stabledt_(
#if (NDIM>2)
             d_omega,
#endif
             dx,xlo,
             ifirst(0),ilast(0),
#if (NDIM>1)
             ifirst(1),ilast(1),
#endif
#if (NDIM>2)
             ifirst(2),ilast(2),
#endif
             ghost_cells(0),
#if (NDIM>1)
             ghost_cells(1),
#endif
#if (NDIM>2)
             ghost_cells(2),
#endif
             uval_current->getPointer(),
             stabdt);

   return stabdt;

}

/********************************************************************
*  计算时间步长.
*********************************************************************/
void RotAdv::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="FLUX" || intc_name=="CONSERVATIVE_DIFFERENCE" || "TAGC");
   #endif

   if(intc_name=="FLUX") {
       computeFluxOnPatch(patch,time,dt);
   } else if(intc_name=="CONSERVATIVE_DIFFERENCE") {
       conservativeDifferenceOnPatch(patch,time,dt);
   } else if(intc_name=="TAGC") {
       int tag_index = this->getTagIndex();
       tagCellsOnPatch(patch,time,tag_index,initial_time);
   } else{
           TBOX_ERROR("\n::computeOnPatch() : component is not matched. "<<endl);
   }

}

/********************************************************************
*  初始化数据片.
*********************************************************************/
void RotAdv::initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT");
   #endif

   (void) time;

   if (initial_time) {
       const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
                                    patch.getPatchData(d_uval_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!uval_current.isNull());
#endif
      hier::IntVector<NDIM> ghost_cells = uval_current->getGhostCellWidth();

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      int model_problem=-1;
      if(d_model_problem=="CIRCLE_ROTATION_DELTA") {
                            model_problem=CIRCLE_ROTATION_DELTA;
      } else if(d_model_problem=="CIRCLE_ROTATION_PARAB") {
                            model_problem=CIRCLE_ROTATION_PARAB;
      } else if(d_model_problem=="SQUARE_ROTATION_DELTA") {
                            model_problem=SQUARE_ROTATION_DELTA;
      } else {
        TBOX_ERROR("RotAdv::initializeDataOnPatch() no model input." << endl);
      }

      initrotation_(model_problem,dx,xlo,d_center,d_radius,
                    ifirst(0),ilast(0),
#if (NDIM>1)
                    ifirst(1),ilast(1),
#endif
#if (NDIM>2)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
#if (NDIM>1)
                    ghost_cells(1),
#endif
#if (NDIM>2)
                    ghost_cells(2),
#endif
                    uval_current->getPointer());
   }

}

/********************************************************************
*  填充物理边界条件.
*********************************************************************/
void RotAdv::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TAGC" || intc_name=="FLUX" || intc_name=="INIT");
   #endif

   NULL_USE(fill_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > uval;
   if(intc_name=="FLUX") uval= patch.getPatchData(d_uval_new_id);
   if(intc_name=="TAGC" || intc_name=="INIT") 
                         uval= patch.getPatchData(d_uval_scratch_id);

   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(!uval.isNull());
   #endif
   hier::IntVector<NDIM> ghost_cells = uval->getGhostCellWidth();

   const tbox::Pointer< hier::PatchGeometry<NDIM> > pgeom =
                                       patch.getPatchGeometry();

   const hier::Box<NDIM>& patch_box = patch.getBox();

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);

         uval->fillAll(0.0,fill_box);
      }
   }
}

/********************************************************************
*  输出数据成员到重启动数据库.
*********************************************************************/
void RotAdv::putToDatabase(tbox::Pointer<tbox::Database> db)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putInteger("d_godunov_order", d_godunov_order);
   db->putString("d_corner_transport", d_corner_transport);
   db->putIntegerArray("d_nghosts", (int*)d_nghosts, NDIM);
   db->putIntegerArray("d_fluxghosts", (int*)d_fluxghosts, NDIM);

   db->putDoubleArray("d_grad_tol", d_grad_tol);
   db->putString("d_model_problem", d_model_problem);
   db->putDoubleArray("d_center", d_center,NDIM);
   db->putDouble("d_radius", d_radius);

#if (NDIM>2)
   db->putDoubleArray("d_omega",d_omega,NDIM);
#endif

}

/********************************************************************
*  打印数据成员.
*********************************************************************/
void RotAdv::printClassData(ostream& os) const
{

   int j;

   os << "\nRotAdv::printClassData..." << endl;
   os << "RotAdv: this = " << (RotAdv*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::UniRectangularGridGeometry<NDIM>*)d_grid_geometry << endl;

   os<< "list (variable, context, patch data index) :\n"
     << "   " << d_uval->getName() <<"  " << d_current->getName()
                                  <<"  " << d_uval_current_id << "\n" 
     << "   " << d_uval->getName() <<"  " << d_new->getName()
                                  <<"  " << d_uval_new_id << "\n" 
     << "   " << d_uval->getName() <<"  " << d_scratch->getName()
                                  <<"  " << d_uval_scratch_id << "\n" 
     << "   " << d_flux->getName() <<"  " << d_new->getName()
                                  <<"  " << d_flux_new_id << "\n" 
     << endl;

   os << "d_model_problem =" << d_model_problem
#if (NDIM>2)
      << " with " << "d_omega = " << d_omega
#endif
      << " with " << "d_center = " << d_center
      << " , d_radius = " << d_radius << endl;

   os << "Parameters for numerical method ..." << endl;
   os << endl;
   os << "   d_godunov_order = " << d_godunov_order << endl;
   os << "   d_corner_transport = " << d_corner_transport << endl;
   os << "   d_nghosts = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   
   os << "   Refinement criteria parameters " << endl;
   os << "       refinement_criteria = Gradient Detector. " << endl;
   for (j = 0; j < d_grad_tol.getSize(); j++) {
      os << "       d_grad_tol[" << j << "] = " << d_grad_tol[j] << endl;
   }
   os << endl;

}

// 以下为私有成员函数.

/********************************************************************
*  从输入文件读取数据成员. 
*********************************************************************/
void RotAdv::getFromInput(tbox::Pointer<tbox::Database> db,
                          bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   /*
    * 注意: 在断点续算时, 当前采用的负载平衡方法应与原计算中的负载平衡方法一致,
    * 例如两者均采用非均匀负载平衡方法, 或均采用均匀负载平衡方法.
    */
   if (db->keyExists("godunov_order")) {
      d_godunov_order = db->getInteger("godunov_order");
      if ( (d_godunov_order != 1) &&
           (d_godunov_order != 2) &&
           (d_godunov_order != 4) ) {
         TBOX_ERROR(d_object_name << ": "
                    << "`godunov_order' in input must be 1, 2, or 4." << endl);
      }
   } else {
      d_godunov_order = db->getIntegerWithDefault("d_godunov_order",
                                                   d_godunov_order);
   }

   if (db->keyExists("corner_transport")) {
      d_corner_transport = db->getString("corner_transport");
      if ( (d_corner_transport != "CORNER_TRANSPORT_1") &&
           (d_corner_transport != "CORNER_TRANSPORT_2") ) {
         TBOX_ERROR(d_object_name << ": "
            << "`corner_transport' in input must be either string"
            << " 'CORNER_TRANSPORT_1' or 'CORNER_TRANSPORT_2'." << endl);
      }
   } else {
      d_corner_transport = db->getStringWithDefault("corner_transport",
                                                     d_corner_transport);
   }

   if (db->keyExists("grad_tol")) {
       d_grad_tol = db->getDoubleArray("grad_tol");
   } else {
       TBOX_ERROR(d_object_name << ": "
              << "No key `grad_tol' found in data. " << endl);
   }

   if (db->keyExists("model_problem")) {
       d_model_problem= db->getString("model_problem");
   } else {
       TBOX_ERROR(d_object_name << ": "
              << "`model_problem' in input must be one of three strings "
              << " `CIRCLE_ROTATION_DELTA' or `CIRCLE_ROTATION_PARAB'"
              << " or `SQUARE_ROTATION_DELTA'. " << endl);
   }

#if (NDIM>2)
   if (db->keyExists("omega")) {
       db->getDoubleArray("omega",d_omega,NDIM);
   } else {
       TBOX_ERROR(d_object_name << ": "
              << " no rotation vector omega found. " << endl);
   }
#endif

   if (db->keyExists("center")) {
       db->getDoubleArray("center",d_center,NDIM);
   } else {
       TBOX_ERROR(d_object_name << ": "
              << " No key `center' found in data." << endl);
   }

   if (db->keyExists("radius")) {
       d_radius = db->getDouble("radius");
   } else {
       TBOX_ERROR(d_object_name << ": "
              << " No key `radius' found in data." << endl);
   }

}

/********************************************************************
*  从输入文件读取数据成员. 
*********************************************************************/
void RotAdv::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = d_nghosts;
   db->getIntegerArray("d_nghosts", tmp_nghosts, NDIM);
   if ( !(d_nghosts == CELLG) ) {
      TBOX_ERROR(d_object_name << ": "
                 << "Key data `d_nghosts' in restart file != CELLG." << endl);
   }
   int* tmp_fluxghosts = d_fluxghosts;
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, NDIM);
   if ( !(d_fluxghosts == FLUXG) ) {
      TBOX_ERROR(d_object_name << ": "
         << "Key data `d_fluxghosts' in restart file != FLUXG." << endl);
   }

   d_grad_tol = db->getDoubleArray("d_grad_tol");
   d_model_problem = db->getString("d_model_problem");
   db->getDoubleArray("d_center",d_center,NDIM);
   d_radius = db->getDouble("d_radius");

#if (NDIM>2)
   db->getDoubleArray("d_omega",d_omega,NDIM);
#endif

}

/*
*************************************************************************
*                                                                       *
* 计算通量.
*                                                                       *
*************************************************************************
*/
#if (NDIM==2)
void RotAdv::computeFluxOnPatch(hier::Patch<NDIM>& patch,
                                const double time,
                                const double dt)
{
   (void) time;

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* xlo = pgeom->getXLower();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval_new_id);
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux_new_id);

   /*
    * 验证由积分对象创建的上下文, 影像区宽度是否正确.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!uval.isNull());
   assert(!flux.isNull());
   assert(uval->getGhostCellWidth() == d_nghosts);
   assert(flux->getGhostCellWidth() == d_fluxghosts);
#endif 

   /*
    * 为网格片数据片分配临时空间.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);

   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  uval->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  flux->getPointer(0),
                  flux->getPointer(1)
                  );
 
   if (d_godunov_order >1) {

      /*
       * 创建临时数组, 用于界面追踪.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::Utilities::imax(Mcells, pbox.numberCells(k));
      }

      // 面心量数组.
      tbox::Array<double>  ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double>  ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double>  ttracrgt(2*FACEG+1+Mcells);

      // 中心量数组.
      tbox::Array<double>  ttcelslp(2*CELLG+Mcells);

      /*
       *  界面追踪, 计算初始的w^L 和 w^R 
       *  输入: w^L, w^R (traced_left/right)
       *  输出: w^L, w^R
      */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    Mcells,xlo,dx,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    Mcells,xlo,dx,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

    }  // if (d_godunov_order > 1) ...
 
    /*
     *  计算网格面上的流.
     *  输入: w^L, w^R (traced_left/right)
     *  输出: F (flux)
     */
   fluxcalculation_(dt,1,0,dx,xlo,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));

    /*
     *  重新计算网格面轨迹, 并进行切向校正.
     *  输入: F (flux)
     *  输出: w^L, w^R (traced_left/right)
     */
   fluxcorrec_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),
               dx,
               uval->getPointer(),
               flux->getPointer(0),
               flux->getPointer(1),
               traced_left.getPointer(0),
               traced_left.getPointer(1),
               traced_right.getPointer(0),
               traced_right.getPointer(1));

   /*
    *  基于更新后的轨迹, 重新计算流.
    *  输入: w^L, w^R (traced_left/right)
    *  输出: Output: F (flux)
    */
   fluxcalculation_(dt,0,0,dx,xlo,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));
}
#endif
#if (NDIM==3)
void RotAdv::computeFluxOnPatch(hier::Patch<NDIM>& patch,
                                const double time,
                                const double dt)
{
   (void) time;
   if (d_corner_transport == "CORNER_TRANSPORT_2") {
      compute3DFluxesWithCornerTransport2(patch, dt);
   } else {
      compute3DFluxesWithCornerTransport1(patch, dt);
   }

}
#endif

/*
*************************************************************************
*                                                                       *
* 用改进的Collella角点输运迎风格式,  
* 近似计算流项.
* 例如: 输入值 corner_transport = CORNER_TRANSPORT_1                *
*                                                                       *
*************************************************************************
*/
#if (NDIM==3)
void RotAdv::compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                                 const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval_new_id);
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux_new_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!uval.isNull());
   assert(!flux.isNull());
   assert(uval->getGhostCellWidth() == d_nghosts);
   assert(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * 为临时的网格数据片分配内存空间.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::FaceData<NDIM,double> temp_traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_traced_right(pbox, 1, d_nghosts);

   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  ifirst(2),ilast(2),
                  uval->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_left.getPointer(2),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  traced_right.getPointer(2),
                  flux->getPointer(0),
                  flux->getPointer(1),
                  flux->getPointer(2));

    /*
    * 如果Godunov格式要求更高的精度, 则使用特征线追踪法进行计算
     */
    if (d_godunov_order >1) {
 
       
       /* 用于计算特征轨迹的临时变量. */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::Utilities::imax(Mcells, pbox.numberCells(k));
      }

      // 面心量数组.
      tbox::Array<double>  ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double>  ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double>  ttracrgt(2*FACEG+1+Mcells);

      // 中心量数组.
      tbox::Array<double>  ttcelslp(2*CELLG+Mcells);

       /*
       *  使用特征轨迹计算w^L和w^R的初值.
        *  输入: w^L, w^R (traced_left/right)
        *  输出: w^L, w^R
        */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),
		    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),
		    ifirst(1),ilast(1),
		    ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());
   }

   /*
    * 基于网格面上的数据, 计算初始流.
    * 输入: w^L, w^R (traced_left/right)
    * 输出: F (通量)
    */
   fluxcalculation_(dt,1,0,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));
   /*
    * 在单元面上使用通量校正更新轨迹状态.
    *  并将结果存储在临时向量(如temp_traced_left/right)中.
    *  输入: F (flux), w^L, w^R (traced_left/right)
    *  输出: temp_traced_left/right
    */
   fluxcorrec2d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 xlo,dx,d_omega,1,
                 uval->getPointer(),
                 flux->getPointer(0),
                 flux->getPointer(1),
                 flux->getPointer(2),
                 traced_left.getPointer(0),
                 traced_left.getPointer(1),
                 traced_left.getPointer(2),
                 traced_right.getPointer(0),
                 traced_right.getPointer(1),
                 traced_right.getPointer(2),
                 temp_traced_left.getPointer(0),
                 temp_traced_left.getPointer(1),
                 temp_traced_left.getPointer(2),
                 temp_traced_right.getPointer(0),
                 temp_traced_right.getPointer(1),
                 temp_traced_right.getPointer(2));

   /*
    *  使用更新的轨迹状态重新计算通量.
    *  输入: temp_traced_left/right
    *  输出: temp_flux
    */
   fluxcalculation_(dt,0,1,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    temp_flux.getPointer(0),
                    temp_flux.getPointer(1),
                    temp_flux.getPointer(2),
                    temp_traced_left.getPointer(0),
                    temp_traced_left.getPointer(1),
                    temp_traced_left.getPointer(2),
                    temp_traced_right.getPointer(0),
                    temp_traced_right.getPointer(1),
                    temp_traced_right.getPointer(2));
   /*
    *  基于其它的切向校正通量差分项, 计算面上的轨迹, 并将结果存储在临时向量中,
    *  (如 temp_traced_left/right). 
    *  输入: F (flux), w^L, w^R (traced_left/right)
    *  输出: temp_traced_left/right
    */
   fluxcorrec2d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 xlo,dx,d_omega,-1,
                 uval->getPointer(),
                 flux->getPointer(0),
                 flux->getPointer(1),
                 flux->getPointer(2),
                 traced_left.getPointer(0),
                 traced_left.getPointer(1),
                 traced_left.getPointer(2),
                 traced_right.getPointer(0),
                 traced_right.getPointer(1),
                 traced_right.getPointer(2),
                 temp_traced_left.getPointer(0),
                 temp_traced_left.getPointer(1),
                 temp_traced_left.getPointer(2),
                 temp_traced_right.getPointer(0),
                 temp_traced_right.getPointer(1),
                 temp_traced_right.getPointer(2));

   /*
    *  计算预测通量值, 并将结果存储在通量向量中.
    * 注意: 将结果直接存储在通量向量中只是为了节省存储量. 
    * 此结果并非最终结果, 在后面还需要进行校正.
    *  输入: temp_traced_left/right
    *  输出: flux
    */
   fluxcalculation_(dt,1,0,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    temp_traced_left.getPointer(0),
                    temp_traced_left.getPointer(1),
                    temp_traced_left.getPointer(2),
                    temp_traced_right.getPointer(0),
                    temp_traced_right.getPointer(1),
                    temp_traced_right.getPointer(2));

   /*
    *  用校正通量的切向差分计算单元面上的轨迹向量,
     *  并将结果存储在 w^L (traced_left) 和w^R (traced_right) 向量中.
    *  输入: temp_flux, flux
    *  输出: w^L, w^R (traced_left/right)
    */
   fluxcorrec3d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                xlo,dx,d_omega,
                uval->getPointer(),
                temp_flux.getPointer(0),
                temp_flux.getPointer(1),
                temp_flux.getPointer(2),
                flux->getPointer(0),
                flux->getPointer(1),
                flux->getPointer(2),
                traced_left.getPointer(0),
                traced_left.getPointer(1),
                traced_left.getPointer(2),
                traced_right.getPointer(0),
                traced_right.getPointer(1),
                traced_right.getPointer(2));
   /*
    *  基于校正后的轨迹, 计算最终的通量.
    *  输入:  w^L, w^R (traced_left/right)
    *  输出:  F (flux)
    */
   fluxcalculation_(dt,0,0,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

}
#endif

/*
*************************************************************************
*                                                                       *
* 使用John Trangenstein提出的三维Collella迎风格式计算通量.        *
* 输入参数: corner_transport = "CORNER_TRANSPORT_2"         *
*                                                                       *
*************************************************************************
*/
#if (NDIM==3)
void RotAdv::compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                                 const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval = 
      patch.getPatchData(d_uval_new_id);
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux = 
      patch.getPatchData(d_flux_new_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!uval.isNull());
   assert(!flux.isNull());
   assert(uval->getGhostCellWidth() == d_nghosts);
   assert(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * 创建临时变量.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::CellData<NDIM,double> third_state(pbox, 1, d_nghosts);

   /*
    * 用中心量变量初始化间断两边的数值解(w^R 和 w^L).
    */
   inittraceflux_(ifirst(0),ilast(0),
		  ifirst(1),ilast(1),
		  ifirst(2),ilast(2),
                  uval->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_left.getPointer(2),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  traced_right.getPointer(2),
                  flux->getPointer(0),
                  flux->getPointer(1),
                  flux->getPointer(2));

   /*
    *  用网格面上当前值计算面上的通量.
    *  输入: w^L, w^R (traced_left/right)
    *  输出: F (flux)
    */
   fluxcalculation_(dt,1,1,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   /*
    * 如果Godunov格式要求更高的精度, 则使用特征线追踪法进行计算
    */
   if (d_godunov_order >1) {

      /*
       * 用于计算特征轨迹的临时变量.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::Utilities::imax(Mcells, pbox.numberCells(k));
      }

      // 临时面心量数组
      tbox::Array<double> ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double> ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double> ttracrgt(2*FACEG+1+Mcells);

      // 临时中心量数组
      tbox::Array<double> ttcelslp(2*CELLG+Mcells);

      /*
       *  使用特征轨迹计算w^L和w^R的初值.
       *  输入: w^L, w^R(traced_left/right)
       *  输出: w^L, w^R
       */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    Mcells,xlo,dx,d_omega,d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

   } //  if (d_godunov_order > 1) ...


   for (int idir = 0; idir < NDIM; idir++) {

      /*
       *  在单元中心近似计算轨迹状态(在第 idir 维方向); 
       *  即: "1/3 state traces"
       *  输入:  F (flux)
       *  输出:  third_state
       */
      onethirdstate_(dt,xlo,dx,d_omega,idir,
                     ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                     uval->getPointer(),
                     flux->getPointer(0),
                     flux->getPointer(1),
                     flux->getPointer(2),
                     third_state.getPointer());
      /*
       *    在除第idir之外的两个方向上使用1/3轨迹状态计算通量.
       *    输入:  third_state
       *    输出:   temp_flux (仅限于两个方向)
       */
      fluxthird_(dt,xlo,dx,d_omega,idir,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 uval->getPointer(),
                 third_state.getPointer(),
                 temp_flux.getPointer(0),
                 temp_flux.getPointer(1),
                 temp_flux.getPointer(2));

      /*
       *    在除第idir维之外的两个方向上使用通量差分校正轨迹状态.
       *    输入:  temp_flux
       *    输出:  w^L, w^R (traced_left/right)
       */
      fluxcorrecjt_(dt,xlo,dx,d_omega,idir,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    temp_flux.getPointer(0),
                    temp_flux.getPointer(1),
                    temp_flux.getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   } // 方向循环...
 
   /*
    *  最后, 使用校正的轨迹状态计算通量.  
    *  输入:  w^L, w^R (traced_left/right)
    *  输出:  F (flux)
    */
   fluxcalculation_(dt,0,0,0,dx,xlo,d_omega,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

}
#endif

/*
*************************************************************************
*                                                                       *
* 在单个网格片上, 对通量进行守恒型差分以更新解.          		*
*                                                                       *
*************************************************************************
*/
void RotAdv::conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                           const double time,
                                           const double dt)
{
   (void) time;
   (void) dt;

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > 
                             patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval_new_id);
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux_new_id); 

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!uval.isNull());
   assert(!flux.isNull());
   assert(uval->getGhostCellWidth() == d_nghosts);
   assert(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   consdiff_(ifirst(0),ilast(0),
#if (NDIM>1)
             ifirst(1),ilast(1),
#endif
#if (NDIM>2)
             ifirst(2),ilast(2),
#endif
             dx,
             flux->getPointer(0),
#if (NDIM>1)
             flux->getPointer(1),
#endif
#if (NDIM>2)
             flux->getPointer(2),
#endif
             uval->getPointer());
}

/********************************************************************
*  标记待细化的网格单元.
*********************************************************************/
void RotAdv::tagCellsOnPatch(hier::Patch<NDIM>& patch,
                             const double  tag_time,
                             const int     tag_indx,
                             const bool    initial_time)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(tag_indx>=0);
   #endif

   NULL_USE(tag_time);
   NULL_USE(initial_time);

   const int error_level_number = patch.getPatchLevelNumber();

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> >
            patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   tbox::Pointer< pdat::CellData<NDIM,int> > tags
                     = patch.getPatchData(tag_indx);

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM, double > > uval =
             patch.getPatchData(d_uval_scratch_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!tags.isNull());
   assert(!uval.isNull());
#endif
   hier::IntVector<NDIM> vghost   = uval->getGhostCellWidth();
   hier::IntVector<NDIM> tagghost = tags->getGhostCellWidth();

   const int size = d_grad_tol.getSize();
   const double tol = ( ( error_level_number < size)
          ? d_grad_tol[error_level_number]
          : d_grad_tol[size-1] );

   tags->fillAll(FALSE);

   detectgrad_(
#if (NDIM==2)
         ifirst(0),ilast(0), ifirst(1),ilast(1),
         vghost(0),tagghost(0),d_nghosts(0),
         vghost(1),tagghost(1),d_nghosts(1),
#endif
#if (NDIM==3)
         ifirst(0),ilast(0), ifirst(1),ilast(1),ifirst(2),ilast(2),
         vghost(0),tagghost(0),d_nghosts(0),
         vghost(1),tagghost(1),d_nghosts(1),
         vghost(2),tagghost(2),d_nghosts(2),
#endif

         dx,tol,TRUE,FALSE,
         uval->getPointer(),
         tags->getPointer());
}

#endif
