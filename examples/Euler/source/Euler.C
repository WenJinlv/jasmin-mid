//
// 文件名:  Euler.C
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 170 $
// 修改  :  $Date: 2007-06-27 08:44:22  $
// 描述  :  Euler方程组：单网格片数值计算子程序.
//




#include "Euler.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

using namespace std;

#include "GridGeometry.h"
#include "UniRectangularPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "FaceIndex.h"
#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/IEEE.h"


//检查物理边界条件的条件编译常数
#define CHECK_BDRY_DATA  (0)

//初始化物理边界条件数据的常数, 可用以错误检查
#define BOGUS_BDRY_DATA   (-9999)

//管理物理边界条件数据的子程序
#if (NDIM == 2) 
#include "UniRectangularBoundaryUtilities2.h"
#endif
#if (NDIM == 3)
#include "UniRectangularBoundaryUtilities3.h"
#endif

//声明Fortran数值子程序的头文件.
#include "EulerFort.h"

//解向量中元素个数
#define NEQU         (NDIM + 2) //NDIM个速度分量+ 压力 + 密度

//变量的影像区宽度
#define CELLG           (4) //格式中中心量变量最大的宽度
#define FACEG           (4) //格式中面心量变量最大的宽度
#define FLUXG           (1) //格式中通量型变量最大的宽度

//物理模型问题的常量
#define PIECEWISE_CONSTANT_X        (10) //黎曼问题（X方向）
#define PIECEWISE_CONSTANT_Y        (11) //黎曼问题（Y方向）
#define PIECEWISE_CONSTANT_Z        (12) //黎曼问题（Z方向）
#define SPHERE                      (40) //球内外爆问题
#define STEP                        (80) //前台阶问题
#define CJ_EXPLORE                  (60) //CJ爆轰问题

//Godunov格式中采用的Riemann解法器
#define APPROX_RIEM_SOLVE   (20) // Colella-Glaz-Riemann解法器 
#define EXACT_RIEM_SOLVE    (21) // 精确Riemann解法器
#define HLLC_RIEM_SOLVE     (22) // HLLC-Riemann解法器

//标记网格单元子程序的常量

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

// 表示重启动文件的版本号的常量
#define EULER_VERSION (3)


/*
*************************************************************************
*                                                                       *
*  Euler类的构造函数, 执行以下操作:                                     *
*  (1) 初始化数据成员                                                   *
*  (2) 创建定义物理模型控制方程中的数值解的变量.                        *
*  (3) 如果进行断点续算, 那么调用getFromRestart().                      *
*  (4) 调用getFromInput(), 从指定数据库中读取数据.                      *
*      该操作可能覆盖从重启动文件中读入的数据.                          *
*                                                                       *
*************************************************************************
*/
Euler::Euler(const string& object_name,
             tbox::Pointer<tbox::Database> input_db,
             tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);


   d_grid_geometry = grid_geom;

   /*
    * 为物理参数指定缺省值.
    */

   d_gamma = 1.4;  // 理想气体的质量热容比
   d_riemann_solve = "APPROX_RIEM_SOLVE";  //描述Riemann数值方法的缺省参数.
   d_godunov_order = 1;
   d_corner_transport = "CORNER_TRANSPORT_1";
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(CELLG == FACEG);
#endif
   d_nghosts = hier::IntVector<NDIM>(CELLG);
   d_fluxghosts = hier::IntVector<NDIM>(FLUXG); 

   // 初始化数据成员.

   d_CJ_eng= 5.31658;  //0.5*d_CJ_velocity*d_CJ_velocity/(d_gamma*d_gamma-1)

   /*
    * 问题类型和初始化数据的缺省值.
    */

   d_radius = tbox::IEEE::getSignalingNaN();
   tbox::IEEE::initializeArrayToSignalingNaN(d_center, NDIM); 
   d_density_inside = tbox::IEEE::getSignalingNaN();
   tbox::IEEE::initializeArrayToSignalingNaN(d_velocity_inside, NDIM); 
   d_pressure_inside = tbox::IEEE::getSignalingNaN();
   d_density_outside = tbox::IEEE::getSignalingNaN();   
   tbox::IEEE::initializeArrayToSignalingNaN(d_velocity_outside, NDIM); 
   d_pressure_outside = tbox::IEEE::getSignalingNaN();   

   d_number_of_intervals = 0;
   d_front_position.resizeArray(0);
   d_interval_density.resizeArray(0);
   d_interval_velocity.resizeArray(0);
   d_interval_pressure.resizeArray(0);

   /*
    * 将物理边界条件设为缺省值. 设为BOGUS_BDRY_DATA用以误差检查.
    */

 #if (NDIM == 2) 
   d_master_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
   }
   
   d_master_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_node_bdry_edge.resizeArray(NUM_2D_NODES);

   for (int ni = 0; ni < NUM_2D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_edge_density.resizeArray(NUM_2D_EDGES);
   d_bdry_edge_velocity.resizeArray(NUM_2D_EDGES*NDIM);
   d_bdry_edge_pressure.resizeArray(NUM_2D_EDGES);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_edge_density);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_edge_velocity);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_edge_pressure);
 #endif  // NDIM == 2

 #if (NDIM == 3)
   d_master_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_scalar_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_vector_bdry_face_conds.resizeArray(NUM_3D_FACES);
   for (int fi = 0; fi < NUM_3D_FACES; fi++) {
      d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_edge_bdry_face.resizeArray(NUM_3D_EDGES);
   for (int ei = 0; ei < NUM_3D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_node_bdry_face.resizeArray(NUM_3D_NODES);

   for (int ni = 0; ni < NUM_3D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_face_density.resizeArray(NUM_3D_FACES);
   d_bdry_face_velocity.resizeArray(NUM_3D_FACES*NDIM);
   d_bdry_face_pressure.resizeArray(NUM_3D_FACES);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_face_density);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_face_velocity);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_face_pressure);
 #endif   // NDIM == 3

   /*
    * 从指定的输入数据库/重启动数据库初始化数据成员.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

  

   /*
    * 用输入数据库/重启动数据库中的数据初始化与问题相关的数据成员.
    */

   if (d_riemann_solve == "APPROX_RIEM_SOLVE") {
      d_riemann_solve_int = APPROX_RIEM_SOLVE;
   } else if (d_riemann_solve == "EXACT_RIEM_SOLVE") {
      d_riemann_solve_int = EXACT_RIEM_SOLVE;
   } else if (d_riemann_solve == "HLLC_RIEM_SOLVE") {
      d_riemann_solve_int = HLLC_RIEM_SOLVE;
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_riemann_solve string = "
         << d_riemann_solve << " encountered in constructor" << endl);
   }

   if (d_data_problem == "PIECEWISE_CONSTANT_X") {
      d_data_problem_int = PIECEWISE_CONSTANT_X;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
      d_data_problem_int = PIECEWISE_CONSTANT_Y;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
      d_data_problem_int = PIECEWISE_CONSTANT_Z;
   } else if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
   } else if (d_data_problem == "STEP") {
      d_data_problem_int = STEP;
   } else if (d_data_problem == "CJ_EXPLORE") {
      d_data_problem_int = CJ_EXPLORE;
      if(d_model_problem=="TWO_POINT_DETONATION")
             d_model_type = 0;
      else if(d_model_problem=="PLATE_DETONATION")
             d_model_type = 1;
      else
           TBOX_ERROR("No model problem:" << d_model_problem <<" input +++ " << endl);
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_data_problem string = "
         << d_data_problem << " encountered in constructor" << endl);
   }



   /*
    * 用输入数据库/重启动数据库的数据后处理物理边界条件	 
    */
 #if (NDIM == 2) 
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == REFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_2D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_edge[i] =
            appu::UniRectangularBoundaryUtilities2::getEdgeLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }

 #endif   // NDIM == 2

 #if (NDIM == 3)
   for (int i = 0; i < NUM_3D_FACES; i++) {
      d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
      d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];

      if (d_master_bdry_face_conds[i] == REFLECT_BC) {
         d_scalar_bdry_face_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_3D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
         d_edge_bdry_face[i] =
            appu::UniRectangularBoundaryUtilities3::getFaceLocationForEdgeBdry(
                                            i, d_master_bdry_edge_conds[i]);
      }
   }

   for (int i = 0; i < NUM_3D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_face[i] =
            appu::UniRectangularBoundaryUtilities3::getFaceLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }

 #endif  //NDIM == 3

 stufprobc_(APPROX_RIEM_SOLVE,EXACT_RIEM_SOLVE,HLLC_RIEM_SOLVE,
           PIECEWISE_CONSTANT_X,PIECEWISE_CONSTANT_Y,PIECEWISE_CONSTANT_Z,
           SPHERE,STEP,CJ_EXPLORE,CELLG,FACEG,FLUXG);

   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();

   d_godunov_flux = 
       new appu::EulerEquationGodunov<NDIM>("Godunov", d_riemann_solve , 
                                 d_gamma, d_godunov_order,d_corner_transport);

}

/*
*************************************************************************
*    Euler 类的空析构函数.	 	                                *
*************************************************************************
*/

Euler::~Euler() 
{ 

tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
} 
/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/

void Euler::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   d_density     = new pdat::CellVariable<NDIM,double>("density", 1);
   d_velocity    = new pdat::CellVariable<NDIM,double>("velocity", NDIM);
   d_pressure    = new pdat::CellVariable<NDIM,double>("pressure", 1);
   d_flux        = new pdat::FaceVariable<NDIM,double>("flux", NEQU);
   d_df          = new pdat::CellVariable<NDIM,double>("CJ_function", 1);

   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");
   d_plot_context = d_current;

   const hier::IntVector<NDIM> zeroghosts(0),oneghosts(1);

   // 存储当前时刻的变量值: (density,current), 影像区宽度为0.
   d_density_current_id = variable_db->registerVariableAndContext(d_density,
                                                           d_current,
                                                           zeroghosts);

   // 存储新时刻的变量值: (density,new), 影像区宽度为格式宽度.
   d_density_new_id = variable_db->registerVariableAndContext(d_density,
                                                           d_new,
                                                           d_nghosts);
   // 存储变量的演算值: (density,scratch), 影像区宽度为等于1.
   d_density_scratch_id = variable_db->registerVariableAndContext(d_density,
                                                           d_scratch,
                                                           oneghosts);

   // 存储当前时刻的变量值: (velocity,current), 影像区宽度为0.
  d_velocity_current_id = 
                         variable_db->registerVariableAndContext(d_velocity,
                                                           d_current,
                                                           zeroghosts);

   // 存储新时刻的变量值: (velocity,new), 影像区宽度为格式宽度.
   d_velocity_new_id = variable_db->registerVariableAndContext(d_velocity,
                                                           d_new,
                                                           d_nghosts);
   // 存储变量的演算值: (velocity,scratch), 影像区宽度为等于1.
   d_velocity_scratch_id = 
                        variable_db->registerVariableAndContext(d_velocity,
                                                           d_scratch,
                                                           oneghosts);

   // 存储当前时刻的变量值: (pressure ,current), 影像区宽度为0.
   d_pressure_current_id = 
                       variable_db->registerVariableAndContext(d_pressure,
                                                           d_current,
                                                           zeroghosts);

   // 存储新时刻的变量值: (pressure,new), 影像区宽度为格式宽度.
   d_pressure_new_id = variable_db->registerVariableAndContext(d_pressure,
                                                           d_new,
                                                           d_nghosts);
   // 存储变量的演算值: (pressure,scratch), 影像区宽度为等于1.
   d_pressure_scratch_id = 
                       variable_db->registerVariableAndContext(d_pressure,
                                                           d_scratch,
                                                           oneghosts);
   // 存储新时刻的通量值: (flux,new), 影像区宽度等于1.
   d_flux_new_id = variable_db->registerVariableAndContext(d_flux,
                                                           d_new,
                                                           d_fluxghosts);
                                                           
   // 存储当前时刻的变量值: (CJ_function ,current), 影像区宽度为0.
   d_df_current_id = variable_db->registerVariableAndContext(d_df ,
                                                           d_current,
                                                           zeroghosts);

   // 存储新时刻的变量值: (CJ_function,new), 影像区宽度为格式宽度.
   d_df_new_id = variable_db->registerVariableAndContext(d_df,
                                                           d_new,
                                                           d_nghosts);
   // 存储变量的演算值: (CJ_function,scratch), 影像区宽度为等于1.
   d_df_scratch_id = variable_db->registerVariableAndContext(d_df,
                                                           d_scratch,
                                                           oneghosts);
        

}
/*
*************************************************************************
*                                                                       *
* 注册JaVis数据输出器, 以采用JaVis工具对绘图文件进行后处理.           *
*                                                                       *
*************************************************************************
*/
void Euler::registerJaVisDataWriter(
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(javis_writer.isNull()));
#endif
   d_javis_writer = javis_writer;

   // 注册可视化量.
   #ifdef HAVE_HDF5
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();
   if (!(d_javis_writer.isNull())) {
      d_javis_writer->registerPlotQuantity("Density",
                                           "SCALAR",
                                           variable_db->mapVariableAndContextToIndex(
                                           d_density, d_plot_context));
      
      d_javis_writer->registerPlotQuantity("Velocity",
                                           "VECTOR",
                                           variable_db->mapVariableAndContextToIndex(
                                           d_velocity, d_plot_context));

      d_javis_writer->registerPlotQuantity("Pressure",
                                           "SCALAR",
                                           variable_db->mapVariableAndContextToIndex(
                                           d_pressure, d_plot_context));
      
      d_javis_writer->registerDerivedPlotQuantity("Total Energy", 
                                                  "SCALAR",
                                                  this);
      d_javis_writer->registerDerivedPlotQuantity("Momentum", 
                                                  "VECTOR",
                                                  this);

                                           
      if (d_data_problem == "CJ_EXPLORE"){
      d_javis_writer->registerPlotQuantity("CJ_func",
                                           "SCALAR",
   variable_db->mapVariableAndContextToIndex(d_df, d_plot_context));
      }
   }
   #endif

}
	 
/********************************************************************
*  初始化积分构件.
*********************************************************************/
void Euler::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件.
        intc->registerInitPatchData(d_density_current_id,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                    d_density_scratch_id);
        intc->registerInitPatchData(d_velocity_current_id,
                                    "USER_DEFINED_REFINE",
                                    d_velocity_scratch_id);
        intc->registerInitPatchData(d_pressure_current_id,
                                    "USER_DEFINED_REFINE",
                                    d_pressure_scratch_id);  
                                      
   }else if(intc_name=="FLUX") { // 计算通量的数值构件.
        intc->registerRefinePatchData(d_density_new_id,
                                      d_density_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE",
                                      d_density_new_id,
                                      d_density_current_id,
                                      d_density_new_id);
        intc->registerRefinePatchData(d_velocity_new_id,
                                      d_velocity_current_id,
                                      "USER_DEFINED_REFINE",
                                      d_velocity_new_id,
                                      d_velocity_current_id,
                                      d_velocity_new_id);
        intc->registerRefinePatchData(d_pressure_new_id,
                                      d_pressure_current_id,
                                      "USER_DEFINED_REFINE",
                                      d_pressure_new_id,
                                      d_pressure_current_id,
                                      d_pressure_new_id);      
                                      
   }else if(intc_name=="SYNC_FLUX") {//沿粗细网格层交接面校正通量的同步构件.
        intc->registerCoarsenPatchData(d_flux_new_id,
                                       "CONSERVATIVE_COARSEN");
   }else if(intc_name=="SYNC_OVERLAY") { // 在细网格层覆盖的内部区域, 同步.
        intc->registerCoarsenPatchData(d_density_new_id,
                                       "CONSERVATIVE_COARSEN",
                                       d_density_new_id);
        intc->registerCoarsenPatchData(d_velocity_new_id,
                                      "USER_DEFINED_COARSEN",
                                       d_velocity_new_id);
        intc->registerCoarsenPatchData(d_pressure_new_id,
                                       "USER_DEFINED_COARSEN",
                                       d_pressure_new_id); 
                                       
   }else if(intc_name=="TAGC"){  // 标记构件.
        intc->registerRefinePatchData(d_density_scratch_id,
                                      d_density_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE");
        intc->registerRefinePatchData(d_velocity_scratch_id,
                                      d_velocity_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE");
        intc->registerRefinePatchData(d_pressure_scratch_id,
                                      d_pressure_current_id,
                                       "CONSERVATIVE_LINEAR_REFINE");  

   }else if(intc_name=="RESET_SOLUTION") { // 接收数值解.
        intc->registerCopyPatchData(d_density_current_id,
                                    d_density_new_id);
        intc->registerCopyPatchData(d_velocity_current_id,
                                    d_velocity_new_id);
        intc->registerCopyPatchData(d_pressure_current_id,
                                    d_pressure_new_id);  
                                    
   }else if(intc_name=="SYNC_COPY") { // 同步前复制数据片
        intc->registerCopyPatchData(d_density_new_id,
                                    d_density_current_id);
        intc->registerCopyPatchData(d_velocity_new_id,
                                    d_velocity_current_id);
        intc->registerCopyPatchData(d_pressure_new_id,
                                    d_pressure_current_id);  
                                    
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 为新值数据片调度内存空间.
        intc->registerPatchData(d_density_new_id);
        intc->registerPatchData(d_velocity_new_id);
        intc->registerPatchData(d_pressure_new_id);
        intc->registerPatchData(d_flux_new_id);    
        
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 为演算数据片调度内存空间.
        intc->registerPatchData(d_density_scratch_id);
        intc->registerPatchData(d_velocity_scratch_id);
        intc->registerPatchData(d_pressure_scratch_id);    
        
   }else if(intc_name!= "TIME_STEP_SIZE" &&
            intc_name!="CONSERVATIVE_DIFFERENCE" ){
        TBOX_ERROR("\n::initializeComponent() : component "
                   << intc_name <<" is not matched. "<<endl);
   }    

 if (d_data_problem == "CJ_EXPLORE"){

   if(intc_name=="INIT") {  // 初值构件.
        intc->registerInitPatchData(d_df_current_id,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                    d_df_scratch_id);
                                      
   }else if(intc_name=="FLUX") { // 计算通量的数值构件.
        intc->registerRefinePatchData(d_df_new_id,
                                      d_df_current_id,
                                      "CONSERVATIVE_LINEAR_REFINE",
                                      d_df_new_id,
                                      d_df_current_id,
                                      d_df_new_id);
   }else if(intc_name=="SYNC_OVERLAY") { // 在细网格层覆盖的内部区域, 同步.
        intc->registerCoarsenPatchData(d_df_new_id,
                                       "CONSERVATIVE_COARSEN",
                                       d_df_new_id);
   }else if(intc_name=="TAGC"){  // 标记构件.
        intc->registerRefinePatchData(d_df_scratch_id,
                                      d_df_current_id,
                                       "CONSERVATIVE_LINEAR_REFINE"); 
   }else if(intc_name=="RESET_SOLUTION") { // 接收数值解. 
        intc->registerCopyPatchData(d_df_current_id,
                                    d_df_new_id);
   }else if(intc_name=="SYNC_COPY") { // 同步前复制数据片
        intc->registerCopyPatchData(d_df_new_id,
                                    d_df_current_id);
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 为新值数据片调度内存空间.
        intc->registerPatchData(d_df_new_id);
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){//为演算数据片调内存空间. 
        intc->registerPatchData(d_df_scratch_id);
   }

 }


}


/*
*************************************************************************
* 初始化网格片内部的数据片.                                             *           
*************************************************************************
*/
void Euler::initializePatchData( hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const bool    initial_time,
                                 const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT");
   #endif

   (void) time;


 if (initial_time) {

      const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const double* xhi = pgeom->getXUpper();

      tbox::Pointer< pdat::CellData<NDIM,double> > density_current = 
          patch.getPatchData(d_density_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > velocity_current = 
          patch.getPatchData(d_velocity_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > pressure_current = 
          patch.getPatchData(d_pressure_current_id);       

      tbox::Pointer< pdat::CellData<NDIM,double> > df_current ; 
      if (d_data_problem == "CJ_EXPLORE")
          df_current = patch.getPatchData(d_df_current_id);


#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!density_current.isNull());
      assert(!velocity_current.isNull());
      assert(!pressure_current.isNull());
#endif
      hier::IntVector<NDIM> ghost_cells = density_current->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(velocity_current->getGhostCellWidth() == ghost_cells);
      assert(pressure_current->getGhostCellWidth() == ghost_cells);
#endif

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();


      if (d_data_problem == "SPHERE") {

         eulerinitsphere_(d_data_problem_int,
                          dx,xlo,xhi,
                          ifirst(0),ilast(0),
                          ifirst(1),ilast(1),
#if (NDIM>2)
                          ifirst(2),ilast(2),
#endif
                          ghost_cells(0),
                          ghost_cells(1),
#if (NDIM>2)
                          ghost_cells(2),
#endif
                          d_gamma,
                          density_current->getPointer(),
                          velocity_current->getPointer(),
                          pressure_current->getPointer(),
                          d_density_inside,
                          d_velocity_inside,
                          d_pressure_inside,
                          d_density_outside,
                          d_velocity_outside,
                          d_pressure_outside,
                          d_center, d_radius);

      } else if (d_data_problem == "CJ_EXPLORE"){

       cjeulerinit_(d_model_type,xlo,dx,
                   ifirst(0),ilast(0),
                   ifirst(1),ilast(1),
#if (NDIM>2)
                   ifirst(2),ilast(2),
#endif
                   ghost_cells(0),
                   ghost_cells(1),
#if (NDIM>2)
                   ghost_cells(2),
#endif
                   density_current->getPointer(),
                   velocity_current->getPointer(),
                   pressure_current->getPointer(),
                   df_current->getPointer(),
                   d_left_point,d_right_point,
                   d_DJ_density,
                   d_CJ_velocity,
                   d_CJ_density,
                   d_CJ_pressure);


       } else{   

         eulerinit_(d_data_problem_int,
                    dx,xlo,xhi,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM>2)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM>2)
                    ghost_cells(2),
#endif
                    d_gamma,
                    density_current->getPointer(),
                    velocity_current->getPointer(),
                    pressure_current->getPointer(),
                    d_number_of_intervals,
                    d_front_position.getPointer(),
                    d_interval_density.getPointer(),
                    d_interval_velocity.getPointer(),
                    d_interval_pressure.getPointer());

      }

  }



}

/*
*************************************************************************
* 计算并返回网格片上的稳定时间步长.                                     *
*************************************************************************
*/

double Euler::getPatchDt(hier::Patch<NDIM>& patch,
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



   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
 
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   const tbox::Pointer< pdat::CellData<NDIM,double> > density_current  = 
      patch.getPatchData(d_density_current_id);
   const tbox::Pointer< pdat::CellData<NDIM,double> > velocity_current = 
      patch.getPatchData(d_velocity_current_id);
   const tbox::Pointer< pdat::CellData<NDIM,double> > pressure_current = 
      patch.getPatchData(d_pressure_current_id); 

   tbox::Pointer< pdat::CellData<NDIM,double> > df_current ; 
      if (d_data_problem == "CJ_EXPLORE")
          df_current = patch.getPatchData(d_df_current_id);
/*
   if (d_data_problem == "CJ_EXPLORE")
   const tbox::Pointer< pdat::CellData<NDIM,double> > df_current = 
      patch.getPatchData(d_df_current_id);
*/
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!density_current.isNull());
   assert(!velocity_current.isNull());
   assert(!pressure_current.isNull());
#endif

   hier::IntVector<NDIM> ghost_cells = density_current->getGhostCellWidth();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(velocity_current->getGhostCellWidth() == ghost_cells);
   assert(pressure_current->getGhostCellWidth() == ghost_cells);
#endif 


   double stabdt = 0.;
   if (d_data_problem == "CJ_EXPLORE"){
   cjstabledt_(dx,
             ifirst(0),ilast(0),
             ifirst(1),ilast(1),
#if (NDIM>2)
             ifirst(2),ilast(2),
#endif
             ghost_cells(0),
             ghost_cells(1),
#if (NDIM>2)
             ghost_cells(2),
#endif
             d_gamma,d_CJ_velocity,
             density_current->getPointer(),
             velocity_current->getPointer(),
             pressure_current->getPointer(),
             df_current->getPointer(),
             stabdt);

     }else{

     stabdt = d_godunov_flux->getGodunovDt(patch, d_density_current_id,
                    d_velocity_current_id,
                    d_pressure_current_id,
                    time);
     }

     return stabdt;

}

/********************************************************************
*  数值计算.
*********************************************************************/
void Euler::computeOnPatch(hier::Patch<NDIM>& patch,
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
          if (d_data_problem == "CJ_EXPLORE"){
          CJ_tagCellsOnPatch(patch,time,tag_index,initial_time);
          }else{
             tagCellsOnPatch(patch,time,tag_index,initial_time);
          }
   } else{
           TBOX_ERROR("\n::computePatch() : component is not matched. "<<endl);
   }

}

/*
*************************************************************************
* 计算网格片上的通量.           					*
*************************************************************************
*/

void Euler::computeFluxOnPatch(hier::Patch<NDIM>& patch,
                                 const double time, 
                                 const double dt)
{
   (void) time;

      d_godunov_flux->getGodunovFluxOnPatch(patch, d_density_new_id,
                                    d_velocity_new_id,
                                    d_pressure_new_id,
                                    d_flux_new_id,
                                    dt);

}


/*
*************************************************************************
* 在单个网格片上, 对通量进行守恒型差分以更新解.                         *
************************************************************************
*/

void Euler::conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                          const double time,
                                          const double dt)
{
   (void) time;
   (void) dt;


   const tbox::Pointer< geom::UniRectangularPatchGeometry<NDIM> > patch_geom = 
   patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();



   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=pbox.lower();
   const hier::Index<NDIM> ilast =pbox.upper();
 
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux        = 
      patch.getPatchData(d_flux_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > density     = 
      patch.getPatchData(d_density_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity    = 
      patch.getPatchData(d_velocity_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure    = 
      patch.getPatchData(d_pressure_new_id);     

   tbox::Pointer< pdat::CellData<NDIM,double> > df ;      
   if (d_data_problem == "CJ_EXPLORE")
      df =  patch.getPatchData(d_df_new_id);
      
   double CJ_dt;


#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!density.isNull());
   assert(!velocity.isNull());
   assert(!pressure.isNull());
   assert(!flux.isNull());
   assert(density->getGhostCellWidth() == d_nghosts);
   assert(velocity->getGhostCellWidth() == d_nghosts);
   assert(pressure->getGhostCellWidth() == d_nghosts);
   assert(flux->getGhostCellWidth() == d_fluxghosts);
#endif

    if (d_data_problem == "CJ_EXPLORE"){
     
#if (NDIM==2) 
    CJ_dt = (d_CJ_rb/d_CJ_velocity)*
                   dx[0]*dx[1]/sqrt(dx[0]*dx[0]+dx[1]*dx[1]);

    cjconsdiff_(d_model_type,xlo,dx,
             ifirst(0),ilast(0),ifirst(1),ilast(1),
             flux->getPointer(0),
             flux->getPointer(1),
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer(),
             df->getPointer(),
             d_gamma,d_DJ_density,d_CJ_density,d_CJ_velocity,
             d_CJ_eng,CJ_dt,d_CJ_nb,d_left_point,d_right_point,time+dt);


#endif

#if (NDIM==3)
    CJ_dt = (d_CJ_rb/d_CJ_velocity)*dx[0]*dx[1]*dx[2]/
                   sqrt(dx[0]*dx[0]*dx[1]*dx[1]+
                        dx[1]*dx[1]*dx[2]*dx[2]+
                        dx[2]*dx[2]*dx[0]*dx[0]);
    cjconsdiff_(d_model_type,xlo,dx,
             ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
             flux->getPointer(0),
             flux->getPointer(1),
             flux->getPointer(2),
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer(),
             df->getPointer(),
             d_gamma,d_DJ_density,d_CJ_density,d_CJ_velocity,
             d_CJ_eng,CJ_dt,d_CJ_nb,d_left_point,d_right_point,time+dt);
#endif

    }else{

#if (NDIM==2)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             d_gamma,
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer());
#endif  //NDIM == 2
#if (NDIM==3)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             flux->getPointer(2),
             d_gamma,
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer());
#endif // NDIM == 3
    }

}

/*
*************************************************************************
* 在一些特殊情况(比如: 使用反射边界条件)下, 重新设置物理边界条件.	*
*************************************************************************
*/

void Euler::boundaryReset(hier::Patch<NDIM>& patch,
                          pdat::FaceData<NDIM,double>& traced_left,  
                          pdat::FaceData<NDIM,double>& traced_right) const
{
   const hier::Index<NDIM> ifirst  =patch.getBox().lower();
   const hier::Index<NDIM> ilast   =patch.getBox().upper();
   int i,idir;
   bool bdry_cell = true;

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   int num_domain_boxes = domain_boxes.getNumberOfBoxes();
   const double* dx = patch_geom->getDx();
   const double* xpatchhi = patch_geom->getXUpper();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   pdat::CellIndex<NDIM> icell = ifirst;
   hier::BoxArray<NDIM> bdrybox(2*NDIM);
   hier::Index<NDIM> ibfirst = ifirst;
   hier::Index<NDIM> iblast  = ilast;
   int bdry_case;

   for (idir=0;idir<NDIM; idir++) {
      ibfirst(idir) = ifirst(idir) -1;
      iblast(idir)  = ifirst(idir) -1;
      bdrybox.getBox(2*idir) = hier::Box<NDIM>(ibfirst,iblast);

      ibfirst(idir) = ilast(idir) +1;
      iblast(idir)  = ilast(idir) +1;
      bdrybox.getBox(2*idir+1) = hier::Box<NDIM>(ibfirst,iblast);
   }

   for (idir=0;idir<NDIM; idir++) {
      int bside = 2*idir;
#if (NDIM == 2) 
      bdry_case = d_master_bdry_edge_conds[bside];
#endif  //NDIM == 2
#if (NDIM == 3)
      bdry_case = d_master_bdry_face_conds[bside];
#endif  //NDIM == 3
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bside)); ic; ic++) {
            for ( i = 0; i < num_domain_boxes; i++ ) {
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
              pdat::FaceIndex<NDIM> sidein = pdat::FaceIndex<NDIM>(ic(),idir,1);
              (traced_left)(sidein,0) = (traced_right)(sidein,0);
            }
         }
      }

      int bnode = 2*idir+1;
#if (NDIM == 2) 
      bdry_case = d_master_bdry_edge_conds[bnode];
#endif
#if (NDIM == 3)
      bdry_case = d_master_bdry_face_conds[bnode];
#endif
// BEGIN SIMPLE-MINDED FIX FOR STEP PROBLEM
      if ( (d_data_problem == "STEP") && (bnode == 1)
          && (fabs(xpatchhi[0]-xdomainhi[0]) < dx[0]) ) {
          bdry_case = FLOW_BC;
      }
// END SIMPLE-MINDED FIX FOR STEP PROBLEM
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bnode)); ic; ic++) {
            for ( i = 0; i < num_domain_boxes; i++ ) {
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
              pdat::FaceIndex<NDIM> sidein = pdat::FaceIndex<NDIM>(ic(),idir,0);
              (traced_right)(sidein,0) = (traced_left)(sidein,0);
            }
         }
      }
   }
}

hier::IntVector<NDIM> Euler::getRefineOpStencilWidth() const
{
   return(hier::IntVector<NDIM>(1));
}

void Euler::postprocessRefineOperator(hier::Patch<NDIM>& fine,
                                  const hier::Patch<NDIM>& coarse,
                                  const hier::Box<NDIM>& fine_box,
                                  const hier::IntVector<NDIM>& ratio,
                                  const string& intc_name)

{


   if(intc_name=="TAGC") return;

   tbox::Pointer< pdat::CellData<NDIM,double> > cdensity;
   tbox::Pointer< pdat::CellData<NDIM,double> > fdensity;
   if(intc_name=="FLUX") {
       cdensity     =  coarse.getPatchData(d_density_new_id);
       fdensity     =  fine.getPatchData(d_density_new_id); 
   }
   if(intc_name=="INIT") {
       cdensity     =  coarse.getPatchData(d_density_scratch_id);
       fdensity     =  fine.getPatchData(d_density_scratch_id);
   }

   tbox::Pointer< pdat::CellData<NDIM,double> > cvelocity;
   tbox::Pointer< pdat::CellData<NDIM,double> > fvelocity;
   if(intc_name=="FLUX") {
      cvelocity     =  coarse.getPatchData(d_velocity_new_id);
      fvelocity     =  fine.getPatchData(d_velocity_new_id);  
   }
   if(intc_name=="INIT") {
      cvelocity     =  coarse.getPatchData(d_velocity_scratch_id);
      fvelocity     = fine.getPatchData(d_velocity_scratch_id); 
   }

   tbox::Pointer< pdat::CellData<NDIM,double> > cpressure;
   tbox::Pointer< pdat::CellData<NDIM,double> > fpressure;
   if(intc_name=="FLUX") {
      cpressure     =  coarse.getPatchData(d_pressure_new_id);
      fpressure     =  fine.getPatchData(d_pressure_new_id);
   }  
   if(intc_name=="INIT") {
      cpressure     =   coarse.getPatchData(d_pressure_scratch_id);
      fpressure     =   fine.getPatchData(d_pressure_scratch_id);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!cdensity.isNull());
   assert(!cvelocity.isNull());
   assert(!cpressure.isNull());
   assert(!fdensity.isNull());
   assert(!fvelocity.isNull());
   assert(!fpressure.isNull());

   hier::IntVector<NDIM> gccheck = cdensity->getGhostCellWidth();
   assert(cvelocity->getGhostCellWidth() == gccheck);
   assert(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   assert(fvelocity->getGhostCellWidth() == gccheck);
   assert(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Box<NDIM> cgbox(cdensity->getGhostBox());

   const hier::Index<NDIM> cilo = cgbox.lower();
   const hier::Index<NDIM> cihi = cgbox.upper();
   const hier::Index<NDIM> filo = fdensity->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdensity->getGhostBox().upper();

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > cgeom = coarse.getPatchGeometry();
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > fgeom = fine.getPatchGeometry();

   const hier::Box<NDIM> coarse_box = hier::Box<NDIM>::coarsen(fine_box, ratio);
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf = fine_box.upper();

   const hier::IntVector<NDIM> cons_ghosts(1);
   pdat::CellData<NDIM,double> conserved(coarse_box, 1, cons_ghosts);

   const hier::IntVector<NDIM> tmp_ghosts(0);

   double* diff0 = new double[coarse_box.numberCells(0)+1];
   pdat::CellData<NDIM,double> slope0(coarse_box, 1, tmp_ghosts);

   double* diff1 = new double[coarse_box.numberCells(1)+1];
   pdat::CellData<NDIM,double> slope1(coarse_box, 1, tmp_ghosts);

#if (NDIM == 3)
   double* diff2 = new double[coarse_box.numberCells(2)+1];
   pdat::CellData<NDIM,double> slope2(coarse_box, 1, tmp_ghosts);
#endif

#if (NDIM == 2)
   pdat::CellData<NDIM,double> flat0(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat1(coarse_box, 1, tmp_ghosts);
   int mc = cihi(0)-cilo(0) + 1;
   mc = tbox::Utilities::imax(mc,cihi(1)-cilo(1) + 1);
   mc = mc + 1;
   double* tflat  = new double[mc];
   double* tflat2 = new double[mc];
   double* tsound = new double[mc];
   double* tdensc = new double[mc];
   double* tpresc = new double[mc];
   double* tvelc  = new double[mc];
   conservlinint2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),  
                    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
                    cilo(0),cilo(1),cihi(0),cihi(1),
                    filo(0),filo(1),fihi(0),fihi(1),
                    ratio,
                    cgeom->getDx(),
                    fgeom->getDx(),
                    d_gamma,
                    cdensity->getPointer(),
                    fdensity->getPointer(),
                    cvelocity->getPointer(),
                    cpressure->getPointer(),
                    fvelocity->getPointer(),                  
                    fpressure->getPointer(),
                    conserved.getPointer(),                  
                    tflat,tflat2,tsound,mc,
                    tdensc,tpresc,tvelc,
                    flat0.getPointer(),
                    flat1.getPointer(),
                    diff0,slope0.getPointer(),
                    diff1,slope1.getPointer());

#endif   //NDIM == 2
#if (NDIM == 3)
   pdat::CellData<NDIM,double> flat0(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat1(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat2(coarse_box, 1, tmp_ghosts);
   int mc = cihi(0)-cilo(0) + 1;
   mc = tbox::Utilities::imax(mc,cihi(1)-cilo(1) + 1);
   mc = tbox::Utilities::imax(mc,cihi(2)-cilo(2) + 1);
   mc = mc + 1;
   double* tflat  = new double[mc];
   double* tflat2 = new double[mc];
   double* tsound = new double[mc];
   double* tdensc = new double[mc];
   double* tpresc = new double[mc];
   double* tvelc  = new double[mc];
   conservlinint3d_(ifirstc(0),ifirstc(1),ifirstc(2),        
                    ilastc(0),ilastc(1),ilastc(2),
                    ifirstf(0),ifirstf(1),ifirstf(2),
                    ilastf(0),ilastf(1),ilastf(2),
                    cilo(0),cilo(1),cilo(2),cihi(0),cihi(1),cihi(2),
                    filo(0),filo(1),filo(2),fihi(0),fihi(1),fihi(2),
                    ratio,
                    cgeom->getDx(),
                    fgeom->getDx(),
                    d_gamma,
                    cdensity->getPointer(),
                    fdensity->getPointer(),
                    cvelocity->getPointer(),
                    cpressure->getPointer(),
                    fvelocity->getPointer(),                  
                    fpressure->getPointer(),
                    conserved.getPointer(),                   
                    tflat,tflat2,tsound,mc,
                    tdensc,tpresc,tvelc,
                    flat0.getPointer(),
                    flat1.getPointer(),
                    flat2.getPointer(),
                    diff0,slope0.getPointer(),
                    diff1,slope1.getPointer(),
                    diff2,slope2.getPointer());
#endif   //NDIM == 3
   delete [] tflat;
   delete [] tflat2;
   delete [] tsound;
   delete [] tdensc;
   delete [] tpresc;
   delete [] tvelc;

   delete [] diff0;
   delete [] diff1;
#if (NDIM == 3)
   delete [] diff2;
#endif

}


/*
*************************************************************************
* 通过守恒型地粗化动量和总能量粗化速度和压力                            *
*************************************************************************
*/

hier::IntVector<NDIM> Euler::getCoarsenOpStencilWidth() const
{
   return(hier::IntVector<NDIM>(0));
}

void Euler::postprocessCoarsenOperator(hier::Patch<NDIM>& coarse,
                                   const hier::Patch<NDIM>& fine,
                                   const hier::Box<NDIM>& coarse_box,
                                   const hier::IntVector<NDIM>& ratio,
                                   const string& intc_name)

{



   if(intc_name=="SYNC_FLUX") return;

   tbox::Pointer< pdat::CellData<NDIM,double> > cdensity;
   tbox::Pointer< pdat::CellData<NDIM,double> > fdensity;
   if(intc_name=="SYNC_OVERLAY") {
       cdensity     =  coarse.getPatchData(d_density_new_id);
       fdensity     =  fine.getPatchData(d_density_new_id); 
   }

   tbox::Pointer< pdat::CellData<NDIM,double> > cvelocity;
   tbox::Pointer< pdat::CellData<NDIM,double> > fvelocity;
   if(intc_name=="SYNC_OVERLAY") {
      cvelocity     =  coarse.getPatchData(d_velocity_new_id);
      fvelocity     =  fine.getPatchData(d_velocity_new_id);  
   }


   tbox::Pointer< pdat::CellData<NDIM,double> > cpressure;
   tbox::Pointer< pdat::CellData<NDIM,double> > fpressure;
   if(intc_name=="SYNC_OVERLAY") {
      cpressure     =  coarse.getPatchData(d_pressure_new_id);
      fpressure     =  fine.getPatchData(d_pressure_new_id);
   }  


#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!cdensity.isNull());
   assert(!cvelocity.isNull());
   assert(!cpressure.isNull());
   assert(!fdensity.isNull());
   assert(!fvelocity.isNull());
   assert(!fpressure.isNull());

   hier::IntVector<NDIM> gccheck = cdensity->getGhostCellWidth();
   assert(cvelocity->getGhostCellWidth() == gccheck);
   assert(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   assert(fvelocity->getGhostCellWidth() == gccheck);
   assert(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Index<NDIM> filo = fdensity->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdensity->getGhostBox().upper();
   const hier::Index<NDIM> cilo = cdensity->getGhostBox().lower();
   const hier::Index<NDIM> cihi = cdensity->getGhostBox().upper();

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > fgeom = fine.getPatchGeometry();
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > cgeom = coarse.getPatchGeometry();

   const hier::Box<NDIM> fine_box = hier::Box<NDIM>::refine(coarse_box, ratio);
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf = fine_box.upper();

   const hier::IntVector<NDIM> cons_ghosts(0);
   pdat::CellData<NDIM,double> conserved(fine_box, 1, cons_ghosts);

#if (NDIM == 2)
   conservavg2d_(ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),  
                 ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
                 filo(0),filo(1),fihi(0),fihi(1),
                 cilo(0),cilo(1),cihi(0),cihi(1),
                 ratio,
                 fgeom->getDx(),
                 cgeom->getDx(),
                 d_gamma,
                 fdensity->getPointer(),
                 cdensity->getPointer(),
                 fvelocity->getPointer(),
                 fpressure->getPointer(),
                 cvelocity->getPointer(),                   
                 cpressure->getPointer(),
                 conserved.getPointer());                  
#endif   //NDIM == 2
#if (NDIM == 3)
   conservavg3d_(ifirstf(0),ifirstf(1),ifirstf(2),          
                 ilastf(0),ilastf(1),ilastf(2),
                 ifirstc(0),ifirstc(1),ifirstc(2),
                 ilastc(0),ilastc(1),ilastc(2),
                 filo(0),filo(1),filo(2),fihi(0),fihi(1),fihi(2),
                 cilo(0),cilo(1),cilo(2),cihi(0),cihi(1),cihi(2),
                 ratio,
                 fgeom->getDx(),
                 cgeom->getDx(),
                 d_gamma,
                 fdensity->getPointer(),
                 cdensity->getPointer(),
                 fvelocity->getPointer(),
                 fpressure->getPointer(),
                 cvelocity->getPointer(),                   
                 cpressure->getPointer(),
                 conserved.getPointer());
#endif   //NDIM == 3

}

/*
 ********************************************************************
 * 设置物理边界条件.                                  *
 ********************************************************************
 */ 

void Euler::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TAGC" || intc_name=="FLUX" || intc_name=="INIT");
   #endif

   NULL_USE(fill_time);

   if (d_data_problem == "CJ_EXPLORE"){
           setPhysicalBoundaryConditions_1( 
         patch, fill_time, ghost_width_to_fill,intc_name);
   }else{  
           setPhysicalBoundaryConditions_2( 
         patch, fill_time, ghost_width_to_fill,intc_name);
   }

}

void Euler::setPhysicalBoundaryConditions_1(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TAGC" || intc_name=="FLUX" || intc_name=="INIT");
   #endif

   NULL_USE(fill_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > density;
   if(intc_name=="FLUX") density= patch.getPatchData(d_density_new_id);
   if(intc_name=="TAGC") density= patch.getPatchData(d_density_scratch_id);
   if(intc_name=="INIT") density= patch.getPatchData(d_density_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > velocity;
   if(intc_name=="FLUX") velocity= patch.getPatchData(d_velocity_new_id);
   if(intc_name=="TAGC") velocity= patch.getPatchData(d_velocity_scratch_id);
   if(intc_name=="INIT") velocity= patch.getPatchData(d_velocity_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > pressure;
   if(intc_name=="FLUX") pressure= patch.getPatchData(d_pressure_new_id);
   if(intc_name=="TAGC") pressure= patch.getPatchData(d_pressure_scratch_id);
   if(intc_name=="INIT") pressure= patch.getPatchData(d_pressure_scratch_id);
  
   tbox::Pointer< pdat::CellData<NDIM,double> > df;
   if (d_data_problem == "CJ_EXPLORE"){
       if(intc_name=="FLUX" )
           df= patch.getPatchData(d_df_new_id);
       if(intc_name=="TAGC" || intc_name=="INIT" )
           df= patch.getPatchData(d_df_scratch_id);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!density.isNull());
   assert(!velocity.isNull());
   assert(!pressure.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = density->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(velocity->getGhostCellWidth() == ghost_cells);
   assert(pressure->getGhostCellWidth() == ghost_cells);
#endif


   const tbox::Pointer< hier::PatchGeometry<NDIM> > pgeom =
            patch.getPatchGeometry();

   const hier::Box<NDIM>& patch_box = patch.getBox();

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);

         density->fillAll(d_DJ_density,fill_box);
         velocity->fillAll(0.0,fill_box);
         pressure->fillAll(1.e-32,fill_box);
         df->fillAll(0,fill_box);

      }
   }

}

void Euler::setPhysicalBoundaryConditions_2(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TAGC" || intc_name=="FLUX" || intc_name=="INIT");
   #endif

   NULL_USE(fill_time);

   tbox::Pointer< pdat::CellData<NDIM,double> > density;
   if(intc_name=="FLUX") density= patch.getPatchData(d_density_new_id);
   if(intc_name=="TAGC") density= patch.getPatchData(d_density_scratch_id);
   if(intc_name=="INIT") density= patch.getPatchData(d_density_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > velocity;
   if(intc_name=="FLUX") velocity= patch.getPatchData(d_velocity_new_id);
   if(intc_name=="TAGC") velocity= patch.getPatchData(d_velocity_scratch_id);
   if(intc_name=="INIT") velocity= patch.getPatchData(d_velocity_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > pressure;
   if(intc_name=="FLUX") pressure= patch.getPatchData(d_pressure_new_id);
   if(intc_name=="TAGC") pressure= patch.getPatchData(d_pressure_scratch_id);
   if(intc_name=="INIT") pressure= patch.getPatchData(d_pressure_scratch_id);
   

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!density.isNull());
   assert(!velocity.isNull());
   assert(!pressure.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = density->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(velocity->getGhostCellWidth() == ghost_cells);
   assert(pressure->getGhostCellWidth() == ghost_cells);
#endif

 
#if (NDIM == 2) 

   /*
    * 设置与网格片棱相交的单元的物理边界条件.
    */
   tbox::Array<int> tmp_edge_scalar_bcond(NUM_2D_EDGES);
   tbox::Array<int> tmp_edge_vector_bcond(NUM_2D_EDGES);
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      tmp_edge_scalar_bcond[i] = d_scalar_bdry_edge_conds[i];
      tmp_edge_vector_bcond[i] = d_vector_bdry_edge_conds[i];
   }

   if (d_data_problem == "STEP") {

      const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = 
         patch.getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const double* xpatchhi = patch_geom->getXUpper();
      const double* xdomainhi = d_grid_geometry->getXUpper();

      if (fabs(xpatchhi[0]-xdomainhi[0]) < dx[0]) {
         tmp_edge_scalar_bcond[XHI] = FLOW_BC;
         tmp_edge_vector_bcond[XHI] = FLOW_BC;
      }

   }

  appu::UniRectangularBoundaryUtilities2::
    fillEdgeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_scalar_bcond,
                           d_bdry_edge_density);
  
  appu::UniRectangularBoundaryUtilities2::
    fillEdgeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_vector_bcond,
                           d_bdry_edge_velocity);

  appu::UniRectangularBoundaryUtilities2::
    fillEdgeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                          tmp_edge_scalar_bcond,
                           d_bdry_edge_pressure);

   /*
    *  设置与网格片结点相交的单元的物理边界条件.
    */
  appu::UniRectangularBoundaryUtilities2::
     fillNodeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_density);

  appu::UniRectangularBoundaryUtilities2::
     fillNodeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_node_conds,
                           d_bdry_edge_velocity);

 appu::UniRectangularBoundaryUtilities2::
   fillNodeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_pressure);


#endif // NDIM == 2

#if (NDIM == 3)

   /*
    *  设置与网格片面相交的单元的物理边界条件.
    */
  appu::UniRectangularBoundaryUtilities3::
    fillFaceBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_density);

  appu::UniRectangularBoundaryUtilities3::
    fillFaceBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_face_conds,
                           d_bdry_face_velocity);

  appu::UniRectangularBoundaryUtilities3::
    fillFaceBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_pressure);

   /*
    *  设置与网格片棱相交的单元的物理边界条件.
    */
  appu::UniRectangularBoundaryUtilities3::
    fillEdgeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_density);

  appu::UniRectangularBoundaryUtilities3::
    fillEdgeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_edge_conds,
                           d_bdry_face_velocity);

  appu::UniRectangularBoundaryUtilities3::
    fillEdgeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_pressure);

   /*
    *  设置与网格片结点相交的单元的物理边界条件.
    */
   appu::UniRectangularBoundaryUtilities3::
      fillNodeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                          d_bdry_face_density);

   appu::UniRectangularBoundaryUtilities3::
      fillNodeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_node_conds,
                           d_bdry_face_velocity);

   appu::UniRectangularBoundaryUtilities3::
      fillNodeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_pressure);

#endif // NDIM == 3

}

/*
*************************************************************************
* 使用梯度检测法标记要细化的单元.                                       *
*************************************************************************
*/

void Euler::tagCellsOnPatch(hier::Patch<NDIM>& patch,
                             const double  tag_time,
                             const int     tag_indx,
                             const bool    initial_time)

{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(tag_indx>=0);
   #endif



   const int error_level_number = patch.getPatchLevelNumber();

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   tbox::Pointer< pdat::CellData<NDIM,int> > tags        = patch.getPatchData(tag_indx);

   hier::Box<NDIM> pbox = patch.getBox(); 
   hier::Box<NDIM> pboxm1 = pbox.grow(pbox,-1);
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   /*
    * 构造定界box
    */
   hier::Box<NDIM> domain;
   for( int i=0; i < domain_boxes.getNumberOfBoxes(); i++ ) {
     domain += domain_boxes(i);
   }

   const hier::Index<NDIM> domfirst=domain.lower();
   const hier::Index<NDIM> domlast =domain.upper();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   /*
    * 创建一组临时标记数据片, 并初始化为假值.
    */
 
   tbox::Pointer< pdat::CellData<NDIM,int> > temp_tags  = new pdat::CellData<NDIM,int>(pbox, 1, d_nghosts);
   temp_tags->fillAll(FALSE);

#if (NDIM==2)
   /*
    * 与具体问题相关的细化标准.
    */
   if (d_data_problem == "STEP") { 
      if (error_level_number < 2) { 
         hier::Box<NDIM> tagbox(hier::Index<NDIM>(9,0), hier::Index<NDIM>(9,3));
         if (error_level_number == 1) {
            tagbox.refine(hier::IntVector<NDIM>(2));
         }
         hier::Box<NDIM> ibox = pbox * tagbox;
     
         for (pdat::CellIterator<NDIM> itc(ibox); itc; itc++) {
            (*temp_tags)(itc(),0) = TRUE; 
         }
      }
   }
#endif  //NDIM == 2

   /*
    * 可选的细化标准包括: DENSITY_DEVIATION, DENSITY_GRADIENT
    *                     PRESSURE_DEVIATION, PRESSURE_GRADIENT
    * 遍历所有的细化标准, 根据当前的时间间隔判断所选用细化标准.
    */
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<NDIM, double > > var;     
      int size = 0;
      double tol = 0.;
      double dev = 0.;
      bool time_allowed = false;

      if (ref == "DENSITY_DEVIATION") {
         var = patch.getPatchData(d_density_scratch_id);
         size = d_density_dev_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_dev_tol[error_level_number] 
                 : d_density_dev_tol[size-1] );
         size = d_density_dev.getSize();
         dev = ( ( error_level_number < size) 
                 ? d_density_dev[error_level_number] 
                 : d_density_dev[size-1] );        
         size = d_density_dev_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_dev_time_min[error_level_number] 
                 : d_density_dev_time_min[size-1] );
         size = d_density_dev_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_dev_time_max[error_level_number] 
                 : d_density_dev_time_max[size-1] );
         time_allowed = (time_min <= tag_time) && (time_max > tag_time);
      }

      if (ref == "DENSITY_GRADIENT") {
         var = patch.getPatchData(d_density_scratch_id);
         size = d_density_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_grad_tol[error_level_number] 
                 : d_density_grad_tol[size-1] );
         size = d_density_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_grad_time_min[error_level_number] 
                 : d_density_grad_time_min[size-1] );
         size = d_density_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_grad_time_max[error_level_number] 
                 : d_density_grad_time_max[size-1] );
         time_allowed = (time_min <= tag_time) && (time_max > tag_time);
      }

      if (ref == "PRESSURE_DEVIATION") {
         var = patch.getPatchData(d_pressure_scratch_id);
         size = d_pressure_dev_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_dev_tol[error_level_number] 
                 : d_pressure_dev_tol[size-1] );
         size = d_pressure_dev.getSize();
         dev = ( ( error_level_number < size) 
                 ? d_pressure_dev[error_level_number] 
                 : d_pressure_dev[size-1] );        
         size = d_pressure_dev_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_dev_time_min[error_level_number] 
                 : d_pressure_dev_time_min[size-1] );
         size = d_pressure_dev_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_dev_time_max[error_level_number] 
                 : d_pressure_dev_time_max[size-1] );
         time_allowed = (time_min <= tag_time) && (time_max > tag_time);
      }

      if (ref == "PRESSURE_GRADIENT") {
         var = patch.getPatchData(d_pressure_scratch_id);
         size = d_pressure_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_grad_tol[error_level_number] 
                 : d_pressure_grad_tol[size-1] );
         size = d_pressure_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_grad_time_min[error_level_number] 
                 : d_pressure_grad_time_min[size-1] );
         size = d_pressure_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_grad_time_max[error_level_number] 
                 : d_pressure_grad_time_max[size-1] );
         time_allowed = (time_min <= tag_time) && (time_max > tag_time);
      }

      if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!var.isNull());
#endif

         hier::IntVector<NDIM> vghost = var->getGhostCellWidth();
         hier::IntVector<NDIM> tagghost = tags->getGhostCellWidth();

         if (ref == "DENSITY_DEVIATION" || ref == "PRESSURE_DEVIATION") {

            /*
             * 检查在当前步之前设置的标记. 
             */ 
            for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
               double locden = tol;
            
               if (tbox::Utilities::fabs((*var)(ic())-dev) > locden) {
                  (*temp_tags)(ic(),0) = TRUE; 
               }
            }
         }

         if (ref == "DENSITY_GRADIENT" || ref == "PRESSURE_GRADIENT") {
            detectgrad_(ifirst(0),ilast(0),
                        ifirst(1),ilast(1),
#if (NDIM>2)
                        ifirst(2),ilast(2),
#endif
                        vghost(0),tagghost(0),d_nghosts(0),
                        vghost(1),tagghost(1),d_nghosts(1),
#if (NDIM>2)
                        vghost(2),tagghost(2),d_nghosts(2),
#endif
                        dx,
                        tol,
                        TRUE, FALSE,
                        var->getPointer(),
                        tags->getPointer(),temp_tags->getPointer());
         }

      }  // if time_allowed 

   }  // loop over criteria

   /*
    * 用临时标记更新单元标记.
    */
   for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
      (*tags)(ic(),0) = (*temp_tags)(ic(),0);
   }

}

void Euler::CJ_tagCellsOnPatch(hier::Patch<NDIM>& patch,
                             const double  tag_time,
                             const int     tag_indx,
                             const bool    initial_time)

{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(tag_indx>=0);
   #endif


   const int error_level_number = patch.getPatchLevelNumber();

   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > 
                 patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   tbox::Pointer< pdat::CellData<NDIM, double > > pressure =
             patch.getPatchData(d_pressure_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM, double > > df =
             patch.getPatchData(d_df_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,int> > tags 
                              = patch.getPatchData(tag_indx);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!pressure.isNull());
   assert(!tags.isNull());
   assert(!df.isNull());
#endif

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   hier::IntVector<NDIM> vghost   = pressure->getGhostCellWidth();
   hier::IntVector<NDIM> tagghost = tags->getGhostCellWidth();

   const int size = d_pressure_grad_tol.getSize();
   const double tol = ( ( error_level_number < size)
          ? d_pressure_grad_tol[error_level_number]
          : d_pressure_grad_tol[size-1] );

   tags->fillAll(FALSE);
   
   cjdetectgrad_(
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
         pressure->getPointer(),
         df->getPointer(),
         tags->getPointer());

}
/*
*************************************************************************
* 将网格片上的总能量和动量(Vis绘图导出量)打包到一个双精度缓冲区.        *
*************************************************************************
*/

bool Euler::packDerivedDataIntoDoubleBuffer(
   double* dbuffer,
   const hier::Patch<NDIM>& patch,
   const hier::Box<NDIM>& region,
   const string& variable_name,
   int depth_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((region * patch.getBox()) == region);
#endif

   bool data_on_patch = FALSE;
   
   tbox::Pointer< pdat::CellData<NDIM,double> > density  =
      patch.getPatchData(d_density, d_plot_context);
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity, d_plot_context);
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
      patch.getPatchData(d_pressure, d_plot_context);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!density .isNull());
   assert(!velocity.isNull());
   assert(!pressure.isNull());
   assert(density->getGhostBox() == patch.getBox());
   assert(velocity->getGhostBox() == patch.getBox());
   assert(pressure->getGhostBox() == patch.getBox());
#endif

   const hier::Box<NDIM>& data_box = density->getGhostBox();
   const int box_w0 = region.numberCells(0);
   const int dat_w0 = data_box.numberCells(0);
   const int box_w1 = region.numberCells(1);
#if (NDIM > 2)
   const int dat_w1 = data_box.numberCells(1);
   const int box_w2 = region.numberCells(2);
#endif

   if (variable_name == "Total Energy") {
      const double *const dens = density->getPointer(); 
      const double *const xvel = velocity->getPointer(0); 
      const double *const yvel = velocity->getPointer(1);
#if (NDIM > 2)
      const double *const zvel = velocity->getPointer(2);
#endif
      const double *const pres = pressure->getPointer(); 

      double valinv = 1.0/(d_gamma-1.0);
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

#if (NDIM > 2)
      for (int i2 = 0; i2 < box_w2; i2++) {
#endif
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               double v2norm = pow(xvel[dat_indx], 2.0)
                             + pow(yvel[dat_indx], 2.0)
#if (NDIM > 2)
                             + pow(zvel[dat_indx], 2.0)
#endif
               ;
               double rho = dens[dat_indx];
               double int_energy = 0.0;
               if (rho > 0.0) {
                  int_energy = valinv * pres[dat_indx] / dens[dat_indx];
               }
               dbuffer[buf_b1+i0] = 
                  dens[dat_indx] * (0.5 * v2norm + int_energy);
            } 
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
#if (NDIM > 2)
         dat_b2 += dat_w1 * dat_w0;
      }
#endif

      data_on_patch = TRUE;

   } else if (variable_name == "Momentum") {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(depth_id < NDIM);
#endif		

      const double *const dens = density->getPointer();
      const double *const vel = velocity->getPointer(depth_id);
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

#if (NDIM > 2)
      for (int i2 = 0; i2 < box_w2; i2++) {
#endif
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               dbuffer[buf_b1+i0] = dens[dat_indx] * vel[dat_indx];
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
#if (NDIM > 2)
         dat_b2 += dat_w1 * dat_w0;
      }
#endif

      data_on_patch = TRUE;

   } else {
      TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
             << "\n    unknown variable_name " << variable_name << "\n");
   }

   return(data_on_patch);
   
}

/*
*************************************************************************
* 将Euler类对象的所有数据成员输出到指定的输出流.                        *
*************************************************************************
*/
void Euler::printClassData(ostream &os) const 
{
   int j,k;

   os << "\nEuler::printClassData..." << endl;
   os << "Euler: this = " << (Euler*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::UniRectangularGridGeometry<NDIM>*)d_grid_geometry << endl;

   os << "Parameters for physical problem ..." << endl;
   os << "   d_gamma = " << d_gamma << endl;

   os << "Numerical method description and ghost sizes..." << endl;
   os << "   d_riemann_solve = " << d_riemann_solve << endl;
   os << "   d_riemann_solve_int = " << d_riemann_solve_int << endl;
   os << "   d_godunov_order = " << d_godunov_order << endl;
   os << "   d_corner_transport = " << d_corner_transport << endl;
   os << "   d_nghosts = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;
   os << "   d_data_problem_int = " << d_data_problem_int << endl;
   
   os << "Problem description and initial data..." << endl;
   if (d_data_problem == "CJ_EXPLORE") {
   os << "   d_model_problem = " << d_model_problem << endl;
   os << "       d_left_point = " << d_left_point << endl;
   os << "       d_right_point = " << d_right_point << endl;

   os << "   d_detonation_model  = " << d_detonation_model << endl;
   os << "       d_DJ_density = " << d_DJ_density  << endl;
   os << "       d_CJ_velocity= " << d_CJ_velocity << endl;
   os << "       d_CJ_density = " << d_CJ_density  << endl;
   os << "       d_CJ_pressure= " << d_CJ_pressure << endl;
   }else if(d_data_problem == "SPHERE") {
   os << "       d_radius = " << d_radius << endl;
   os << "       d_center = " ;
      for (j=0;j<NDIM;j++) os << d_center[j]<<" ";
   os << endl;
   os << "       d_density_inside = " << d_density_inside << endl;
   os << "       d_velocity_inside = " ;
     for (j=0;j<NDIM;j++) os << d_velocity_inside[j]<<" ";
   os << endl;
   os << "       d_pressure_inside = " << d_pressure_inside << endl;
   os << "       d_density_outside = " << d_density_outside << endl;
   os << "       d_velocity_outside = " ;
     for (j=0;j<NDIM;j++) os << d_velocity_outside[j]<<" ";
   os << endl;
   os << "       d_pressure_outside = " << d_pressure_outside << endl;
   }else{
   os << "       d_number_of_intervals = " << d_number_of_intervals << endl;
   os << "       d_front_position = ";
   for (k = 0; k < d_number_of_intervals-1; k++) {
      os << d_front_position[k] << "  "; 
   } 
   os << endl;
   os << "       d_interval_density = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            " << d_interval_density[k] << endl;
   } 
   os << "       d_interval_velocity = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            ";
      for (j = 0; j < NDIM; j++) {
         os << d_interval_velocity[k*NDIM+j] << "  "; 
      }
      os << endl;
   } 
   os << "       d_interval_pressure = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            " << d_interval_pressure[k] << endl;
   }
   }

   os << "   Boundary condition data " << endl;

#if (NDIM == 2) 
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      if (d_master_bdry_edge_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_edge_density[" << j << "] = "
            << d_bdry_edge_density[j] << endl;
         os << "         d_bdry_edge_velocity[" << j << "] = "
            << d_bdry_edge_velocity[j*NDIM+0] << " , "
            << d_bdry_edge_velocity[j*NDIM+1] << endl;
         os << "         d_bdry_edge_pressure[" << j << "] = "
            << d_bdry_edge_pressure[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_edge[" << j << "] = "
         << d_node_bdry_edge[j] << endl;
   }
#endif
#if (NDIM == 3)
   for (j = 0; j < d_master_bdry_face_conds.getSize(); j++) {
      os << "\n       d_master_bdry_face_conds[" << j << "] = "
         << d_master_bdry_face_conds[j] << endl;
      os << "       d_scalar_bdry_face_conds[" << j << "] = "
         << d_scalar_bdry_face_conds[j] << endl;
      os << "       d_vector_bdry_face_conds[" << j << "] = "
         << d_vector_bdry_face_conds[j] << endl;
      if (d_master_bdry_face_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_face_density[" << j << "] = "
            << d_bdry_face_density[j] << endl;
         os << "         d_bdry_face_velocity[" << j << "] = "
            << d_bdry_face_velocity[j*NDIM+0] << " , "
            << d_bdry_face_velocity[j*NDIM+1] << " , "
            << d_bdry_face_velocity[j*NDIM+2] << endl;
         os << "         d_bdry_face_pressure[" << j << "] = "
            << d_bdry_face_pressure[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      os << "       d_edge_bdry_face[" << j << "] = "
         << d_edge_bdry_face[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_face[" << j << "] = "
         << d_node_bdry_face[j] << endl;
   }
#endif

   os << "   Refinement criteria parameters " << endl;
 
   for (j = 0; j < d_refinement_criteria.getSize(); j++) {
      os << "       d_refinement_criteria[" << j << "] = "
         << d_refinement_criteria[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_tol.getSize(); j++) {
      os << "       d_density_dev_tol[" << j << "] = "
         << d_density_dev_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev.getSize(); j++) {
      os << "       d_density_dev[" << j << "] = "
         << d_density_dev[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_time_max.getSize(); j++) {
      os << "       d_density_dev_time_max[" << j << "] = "
         << d_density_dev_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_time_min.getSize(); j++) {
      os << "       d_density_dev_time_min[" << j << "] = "
         << d_density_dev_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_tol.getSize(); j++) {
      os << "       d_density_grad_tol[" << j << "] = "
         << d_density_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_time_max.getSize(); j++) {
      os << "       d_density_grad_time_max[" << j << "] = "
         << d_density_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_time_min.getSize(); j++) {
      os << "       d_density_grad_time_min[" << j << "] = "
         << d_density_grad_time_min[j] << endl;
   }
   os << endl;

   for (j = 0; j < d_pressure_dev_tol.getSize(); j++) {
      os << "       d_pressure_dev_tol[" << j << "] = "
         << d_pressure_dev_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev.getSize(); j++) {
      os << "       d_pressure_dev[" << j << "] = "
         << d_pressure_dev[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev_time_max.getSize(); j++) {
      os << "       d_pressure_dev_time_max[" << j << "] = "
         << d_pressure_dev_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev_time_min.getSize(); j++) {
      os << "       d_pressure_dev_time_min[" << j << "] = "
         << d_pressure_dev_time_min[j] << endl;
   }
   os << endl; 
   for (j = 0; j < d_pressure_grad_tol.getSize(); j++) {
      os << "       d_pressure_grad_tol[" << j << "] = "
         << d_pressure_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_grad_time_max.getSize(); j++) {
      os << "       d_pressure_grad_time_max[" << j << "] = "
         << d_pressure_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_grad_time_min.getSize(); j++) {
      os << "       d_pressure_grad_time_min[" << j << "] = "
         << d_pressure_grad_time_min[j] << endl;
   }
   os << endl;

}

/*
*************************************************************************
* 从输入数据库读取数据. 如果从输入数据库和重启动数据库中读取相同数据成  *
* 员的值.那么用输入数据库中的值覆盖重启动数据库的值.                    *
*************************************************************************
*/

void Euler::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("gamma") ) {
      d_gamma = db->getDouble("gamma");
   }

   if (db->keyExists("riemann_solve")) {
      d_riemann_solve = db->getString("riemann_solve");
      if ( (d_riemann_solve != "APPROX_RIEM_SOLVE") &&
           (d_riemann_solve != "EXACT_RIEM_SOLVE") &&
           (d_riemann_solve != "HLLC_RIEM_SOLVE") ) {
         TBOX_ERROR(d_object_name << ": "
                    << "`riemann_solve' in input must be either string "
                    << "'APPROX_RIEM_SOLVE', 'EXACT_RIEM_SOLVE', "
                    << "'HLLC_RIEM_SOLVE'." << endl);

      }
   } else {
      d_riemann_solve = db->getStringWithDefault("d_riemann_solve",
                                                  d_riemann_solve);
   }

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


   if (db->keyExists("Refinement_data")) {
      tbox::Pointer<tbox::Database> refine_db = db->getDatabase("Refinement_data");
      tbox::Array<string> refinement_keys = refine_db->getAllKeys();
      int num_keys = refinement_keys.getSize();

      if (refine_db->keyExists("refine_criteria")) {
         d_refinement_criteria = 
            refine_db->getStringArray("refine_criteria");
      } else {
         TBOX_WARNING(d_object_name << ": "
                   << "No key `refine_criteria' found in data for"
                   << " RefinementData. No refinement will occur." << endl);
      }
       
      tbox::Array<string> ref_keys_defined(num_keys);
      int def_key_cnt = 0;
      tbox::Pointer<tbox::Database> error_db;
      for (int i = 0; i < refinement_keys.getSize(); i++) {
           
         string error_key = refinement_keys[i];
         error_db.setNull();

         if ( !(error_key == "refine_criteria") ) {

            if ( !(error_key == "DENSITY_DEVIATION" ||
                   error_key == "DENSITY_GRADIENT" ||
                   error_key == "PRESSURE_DEVIATION" ||
                   error_key == "PRESSURE_GRADIENT" )){
               TBOX_ERROR(d_object_name << ": "
                         << "Unknown refinement criteria: " << error_key
                         << "\nin input." << endl);
            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               def_key_cnt++;
            }
               
            if (!error_db.isNull() && error_key == "DENSITY_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_density_dev_tol = 
                  error_db->getDoubleArray("dev_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `dev_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("density_dev")) {
                  d_density_dev = 
                  error_db->getDoubleArray("density_dev");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `density_dev' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_dev_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_dev_time_max.resizeArray(1);
                  d_density_dev_time_max[0] = tbox::IEEE::getDBL_MAX();
               }

               if (error_db->keyExists("time_min")){
                  d_density_dev_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_dev_time_min.resizeArray(1);
                  d_density_dev_time_min[0] = 0.;
               } 

            }
              
            if (!error_db.isNull() && error_key == "DENSITY_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_density_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_grad_time_max.resizeArray(1);
                  d_density_grad_time_max[0] = tbox::IEEE::getDBL_MAX();
               }

               if (error_db->keyExists("time_min")){
                  d_density_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_grad_time_min.resizeArray(1);
                  d_density_grad_time_min[0] = 0.;
               } 

            }
              
            if (!error_db.isNull() && error_key == "PRESSURE_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_pressure_dev_tol = 
                  error_db->getDoubleArray("dev_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `dev_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("pressure_dev")) {
                  d_pressure_dev = 
                  error_db->getDoubleArray("pressure_dev");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `pressure_dev' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_dev_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_dev_time_max.resizeArray(1);
                  d_pressure_dev_time_max[0] = tbox::IEEE::getDBL_MAX();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_dev_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_dev_time_min.resizeArray(1);
                  d_pressure_dev_time_min[0] = 0.;
               } 

            }

            if (!error_db.isNull() && error_key == "PRESSURE_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_pressure_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_grad_time_max.resizeArray(1);
                  d_pressure_grad_time_max[0] = tbox::IEEE::getDBL_MAX();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_grad_time_min.resizeArray(1);
                  d_pressure_grad_time_min[0] = 0.;
               }

            }
      
         } 

      } // 遍历所有细化标准     


      /* 
       * 检查所要求的每个细化标准在输入数据库中是否都存在.
       */
      for (int k0 = 0; k0 < d_refinement_criteria.getSize(); k0++) {
         string use_key = d_refinement_criteria[k0];
         bool key_found = false;
         for (int k1 = 0; k1 < def_key_cnt; k1++) {
             string def_key = ref_keys_defined[k1];
             if (def_key == use_key) key_found = true;
         }

         if (!key_found) {
             TBOX_ERROR(d_object_name << ": "
                       << "No input found for specified refine criteria: "
                       << d_refinement_criteria[k0] << endl);
         }
      }

   } // if "Refinement_data" db entry exists

   if (!is_from_restart) { 

      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`data_problem' value not found in input." << endl);
      }

      tbox::Pointer<tbox::Database> init_data_db;
      if (db->keyExists("Initial_data")) {
         init_data_db = db->getDatabase("Initial_data");
      } else {
         TBOX_ERROR(d_object_name << ": "
                  << "No `Initial_data' database found in input." << endl);
      }

      bool found_problem_data = false;

      if (d_data_problem == "SPHERE") {

         if (init_data_db->keyExists("radius")) {
            d_radius = init_data_db->getDouble("radius");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`radius' input required for SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("center")) {
            init_data_db->getDoubleArray("center", d_center, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`center' input required for SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("density_inside")) {
            d_density_inside = init_data_db->getDouble("density_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`density_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("velocity_inside")) {
            init_data_db->getDoubleArray("velocity_inside",
                                         d_velocity_inside, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`velocity_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("pressure_inside")) {
            d_pressure_inside = init_data_db->getDouble("pressure_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`pressure_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("density_outside")) {
            d_density_outside = init_data_db->getDouble("density_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`density_outside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("velocity_outside")) {
            init_data_db->getDoubleArray("velocity_outside",
                                         d_velocity_outside, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`velocity_outside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("pressure_outside")) {
            d_pressure_outside = init_data_db->getDouble("pressure_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`pressure_outside' input required for "
               << "SPHERE problem." << endl);
         }
   
         found_problem_data = true;

      }
      
      if (d_data_problem == "CJ_EXPLORE") {
      
      if (init_data_db->keyExists("model_problem")) {
          d_model_problem = init_data_db->getString("model_problem");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `model_problem' found"<<endl);
      }
      if (init_data_db->keyExists("left_point")) {
          init_data_db->getDoubleArray("left_point",d_left_point,NDIM);
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `left_point' found"<<endl);
      }
      if (init_data_db->keyExists("right_point")) {
          init_data_db->getDoubleArray("right_point",d_right_point,NDIM);
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `right_point' found"<<endl);
      }
      if (init_data_db->keyExists("detonation_model")) {
          d_detonation_model= init_data_db->getString("detonation_model");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `detonation_model' found"<<endl);
      }
      if (init_data_db->keyExists("DJ_density")) {
          d_DJ_density = init_data_db->getDouble("DJ_density");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `DJ_density' found"<<endl);
      }
      if (init_data_db->keyExists("CJ_velocity")) {
          d_CJ_velocity = init_data_db->getDouble("CJ_velocity");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `CJ_velocity' found"<<endl);
      }
      if (init_data_db->keyExists("CJ_density")) {
          d_CJ_density = init_data_db->getDouble("CJ_density");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `CJ_density' found"<<endl);
      }
      if (init_data_db->keyExists("CJ_pressure")) {
          d_CJ_pressure = init_data_db->getDouble("CJ_pressure");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `CJ_pressure' found"<<endl);
      }
      if (init_data_db->keyExists("CJ_rb")) {
          d_CJ_rb= init_data_db->getDouble("CJ_rb");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `CJ_rb' found"<<endl);
      }
      if (init_data_db->keyExists("CJ_nb")) {
          d_CJ_nb= init_data_db->getDouble("CJ_nb");
      } else {
          TBOX_ERROR(d_object_name << ": "
                           << "No key `CJ_nb' found"<<endl);
      }
         found_problem_data = true;

   }   


      if (!found_problem_data &&
          (d_data_problem == "PIECEWISE_CONSTANT_X") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
          (d_data_problem == "STEP")) {

          int idir = 0;
          if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
             idir = 1;
          }

          if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
             if (NDIM < 3) {
                TBOX_ERROR(d_object_name << ": `PIECEWISE_CONSTANT_Z' "
                           << "problem invalid in 2 dimensions." << endl);
             }
             idir = 2;
          }

          tbox::Array<string> init_data_keys = init_data_db->getAllKeys();

          if (init_data_db->keyExists("front_position")) {
             d_front_position = init_data_db->getDoubleArray("front_position");
          } else {
             TBOX_ERROR(d_object_name << ": "
                << "`front_position' input required for "
                << "PIECEWISE_CONSTANT_* problem." << endl);
          }

          d_number_of_intervals = 
             tbox::Utilities::imin(d_front_position.getSize()+1,
                                  init_data_keys.getSize()-1);

          d_front_position.resizeArray(d_front_position.getSize()+1);
          d_front_position[d_front_position.getSize()-1] = 
             d_grid_geometry->getXUpper()[idir]; 

          d_interval_density.resizeArray(d_number_of_intervals);
          d_interval_velocity.resizeArray(d_number_of_intervals*NDIM);
          d_interval_pressure.resizeArray(d_number_of_intervals);

          int i = 0;
          int nkey = 0;
          bool found_interval_data = false;

          while (   !found_interval_data 
                 && (i < d_number_of_intervals)
                 && (nkey < init_data_keys.getSize()) ) {

             if ( !(init_data_keys[nkey] == "front_position") ) {

                tbox::Pointer<tbox::Database> interval_db =
                      init_data_db->getDatabase(init_data_keys[nkey]);

                   readStateDataEntry(interval_db,
                                      init_data_keys[nkey],
                                      i,
                                      d_interval_density,
                                      d_interval_velocity,
                                      d_interval_pressure);

                i++;
   
                found_interval_data = (i == d_number_of_intervals);
   
             }
   
             nkey++;
   
          }

          if (!found_interval_data) {
             TBOX_ERROR(d_object_name << ": "
                        << "Insufficient interval data given in input"
                        << " for PIECEWISE_CONSTANT_* or STEP problem." << endl);
          }
      
          found_problem_data = true;

      }

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
            << "`Initial_data' database found in input." 
            << " But bad data supplied." << endl);
      } 
 
   } // if !is_from_restart read in problem data

   hier::IntVector<NDIM> periodic = d_grid_geometry->getPeriodicShift();
   int num_per_dirs = 0;
   for (int id = 0; id < NDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < NDIM) {

      if (db->keyExists("Boundarydata")&&d_data_problem != "CJ_EXPLORE") {

         tbox::Pointer<tbox::Database> bdry_db = db->getDatabase("Boundarydata");

#if (NDIM == 2)
          appu::UniRectangularBoundaryUtilities2::readBoundaryInput( this,
                                                        bdry_db,
                                                        d_master_bdry_edge_conds,
                                                        d_master_bdry_node_conds,
                                                        periodic);
#endif
#if (NDIM == 3)
         appu::UniRectangularBoundaryUtilities3::readBoundaryInput(this,
                                                        bdry_db,
                                                        d_master_bdry_face_conds,
                                                        d_master_bdry_edge_conds,
                                                        d_master_bdry_node_conds,
                                                        periodic);
#endif

      } else {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `Boundarydata' not found in input. " << endl);
      }

   }

}

/*
*************************************************************************
* 将Euler类的数据成员写入到重启动数据库.                                *
*************************************************************************
*/

void Euler::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putInteger("EULER_VERSION", EULER_VERSION);

   db->putDouble("d_gamma", d_gamma);

   db->putString("d_riemann_solve", d_riemann_solve);
   db->putInteger("d_godunov_order", d_godunov_order);
   db->putString("d_corner_transport", d_corner_transport);
   db->putIntegerArray("d_nghosts", d_nghosts, NDIM);
   db->putIntegerArray("d_fluxghosts", d_fluxghosts, NDIM);

   db->putString("d_data_problem", d_data_problem);
   
   if (d_data_problem == "CJ_EXPLORE") {
   db->putString("d_model_problem", d_model_problem);
   db->putDoubleArray("d_left_point", d_left_point,NDIM);
   db->putDoubleArray("d_right_point", d_right_point,NDIM);
   db->putString("d_detonation_model",d_detonation_model);
   db->putDouble("d_DJ_density", d_DJ_density);
   db->putDouble("d_CJ_velocity", d_CJ_velocity);
   db->putDouble("d_CJ_density", d_CJ_density);
   db->putDouble("d_CJ_pressure", d_CJ_pressure);
   db->putDouble("d_CJ_nb", d_CJ_nb);
   db->putDouble("d_CJ_rb", d_CJ_rb);
   }
   
   if (d_data_problem == "SPHERE") {
      db->putDouble("d_radius", d_radius);
      db->putDoubleArray("d_center", d_center, NDIM);
      db->putDouble("d_density_inside", d_density_inside);
      db->putDoubleArray("d_velocity_inside", d_velocity_inside, NDIM);
      db->putDouble("d_pressure_inside", d_pressure_inside);
      db->putDouble("d_density_outside", d_density_outside);
      db->putDoubleArray("d_velocity_outside", d_velocity_outside, NDIM);
      db->putDouble("d_pressure_outside", d_pressure_outside);
   }

   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "STEP") ) {
      db->putInteger("d_number_of_intervals", d_number_of_intervals);
      if (d_number_of_intervals > 0) {
         db->putDoubleArray("d_front_position", d_front_position);   
         db->putDoubleArray("d_interval_density", d_interval_density);   
         db->putDoubleArray("d_interval_velocity", d_interval_velocity);   
         db->putDoubleArray("d_interval_pressure", d_interval_pressure);   
      }
   }

   db->putIntegerArray("d_master_bdry_edge_conds", d_master_bdry_edge_conds);
   db->putIntegerArray("d_master_bdry_node_conds", d_master_bdry_node_conds);

#if (NDIM == 2) 
   db->putDoubleArray("d_bdry_edge_density", d_bdry_edge_density);
   db->putDoubleArray("d_bdry_edge_velocity", d_bdry_edge_velocity);
   db->putDoubleArray("d_bdry_edge_pressure", d_bdry_edge_pressure);
#endif
#if (NDIM == 3)
   db->putIntegerArray("d_master_bdry_face_conds", d_master_bdry_face_conds);

   db->putDoubleArray("d_bdry_face_density", d_bdry_face_density);
   db->putDoubleArray("d_bdry_face_velocity", d_bdry_face_velocity);
   db->putDoubleArray("d_bdry_face_pressure", d_bdry_face_pressure);
#endif

   if (d_refinement_criteria.getSize() > 0) {
      db->putStringArray("d_refinement_criteria", d_refinement_criteria);
   }
   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         db->putDoubleArray("d_density_dev_tol", 
                            d_density_dev_tol);
         db->putDoubleArray("d_density_dev", 
                            d_density_dev);
         db->putDoubleArray("d_density_dev_time_max", 
                            d_density_dev_time_max);
         db->putDoubleArray("d_density_dev_time_min", 
                            d_density_dev_time_min);

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         db->putDoubleArray("d_density_grad_tol", 
                            d_density_grad_tol);
         db->putDoubleArray("d_density_grad_time_max", 
                            d_density_grad_time_max);
         db->putDoubleArray("d_density_grad_time_min", 
                            d_density_grad_time_min);

      } else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         db->putDoubleArray("d_pressure_dev_tol", 
                            d_pressure_dev_tol);
         db->putDoubleArray("d_pressure_dev", 
                            d_pressure_dev);
         db->putDoubleArray("d_pressure_dev_time_max", 
                            d_pressure_dev_time_max);
         db->putDoubleArray("d_pressure_dev_time_min", 
                            d_pressure_dev_time_min);

      } else if  (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         db->putDoubleArray("d_pressure_grad_tol", 
                            d_pressure_grad_tol);
         db->putDoubleArray("d_pressure_grad_time_max", 
                            d_pressure_grad_time_max);
         db->putDoubleArray("d_pressure_grad_time_min", 
                            d_pressure_grad_time_min);

      } 

   }

}

/*
*************************************************************************
*    从重启动数据库中读取数据.                                          *
*************************************************************************
*/
void Euler::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
   }

   int ver = db->getInteger("EULER_VERSION");
   if (ver != EULER_VERSION) {
      TBOX_ERROR(d_object_name << ": "
          << "Restart file version different than class version." << endl);
   }

   d_gamma = db->getDouble("d_gamma");

   d_riemann_solve = db->getString("d_riemann_solve");
   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = d_nghosts;
   db->getIntegerArray("d_nghosts", tmp_nghosts, NDIM);
   for (int i = 0; i < NDIM; i++) {
      if (d_nghosts(i) != CELLG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_nghosts' in restart file != CELLG." << endl);
      }
   }
   int* tmp_fluxghosts = d_fluxghosts;
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, NDIM);
   for (int i = 0; i < NDIM; i++) {
      if (d_fluxghosts(i) != FLUXG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_fluxghosts' in restart file != FLUXG." << endl);
      }
   }

   d_data_problem = db->getString("d_data_problem");

   if (d_data_problem == "CJ_EXPLORE") {
   d_model_problem = db->getString("d_model_problem");
   db->getDoubleArray("d_left_point", d_left_point,NDIM);
   db->getDoubleArray("d_right_point", d_right_point,NDIM);
   d_detonation_model=db->getString("d_detonation_model");
   d_DJ_density=db->getDouble("d_DJ_density");
   d_CJ_velocity=db->getDouble("d_CJ_velocity");
   d_CJ_density=db->getDouble("d_CJ_density");
   d_CJ_pressure=db->getDouble("d_CJ_pressure");
   d_CJ_nb  =db->getDouble("d_CJ_nb");
   d_CJ_rb  =db->getDouble("d_CJ_rb");
   }
   
   if (d_data_problem == "SPHERE") {
      d_radius = db->getDouble("d_radius");
      db->getDoubleArray("d_center", d_center, NDIM);
      d_density_inside = db->getDouble("d_density_inside");
      db->getDoubleArray("d_velocity_inside", d_velocity_inside, NDIM);
      d_pressure_inside = db->getDouble("d_pressure_inside");
      d_density_outside = db->getDouble("d_density_outside");
      db->getDoubleArray("d_velocity_outside", d_velocity_outside, NDIM);
      d_pressure_outside = db->getDouble("d_pressure_outside");
   } 

   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "STEP") ) {
      d_number_of_intervals = db->getInteger("d_number_of_intervals");
      if (d_number_of_intervals > 0) {
         d_front_position = db->getDoubleArray("d_front_position");
         d_interval_density = db->getDoubleArray("d_interval_density");
         d_interval_velocity = db->getDoubleArray("d_interval_velocity");
         d_interval_pressure = db->getDoubleArray("d_interval_pressure");
      }
   }

   d_master_bdry_edge_conds = db->getIntegerArray("d_master_bdry_edge_conds");
   d_master_bdry_node_conds = db->getIntegerArray("d_master_bdry_node_conds");

#if (NDIM == 2) 
   d_bdry_edge_density = db->getDoubleArray("d_bdry_edge_density");
   d_bdry_edge_velocity = db->getDoubleArray("d_bdry_edge_velocity");
   d_bdry_edge_pressure = db->getDoubleArray("d_bdry_edge_pressure");
#endif
#if (NDIM == 3)
   d_master_bdry_face_conds = db->getIntegerArray("d_master_bdry_face_conds");

   d_bdry_face_density = db->getDoubleArray("d_bdry_face_density");
   d_bdry_face_velocity = db->getDoubleArray("d_bdry_face_velocity");
   d_bdry_face_pressure = db->getDoubleArray("d_bdry_face_pressure");
#endif

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringArray("d_refinement_criteria");
   }

   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         d_density_dev_tol = 
            db->getDoubleArray("d_density_dev_tol");
         d_density_dev = 
            db->getDoubleArray("d_density_dev");
         d_density_dev_time_max = 
            db->getDoubleArray("d_density_dev_time_max");
         d_density_dev_time_min = 
            db->getDoubleArray("d_density_dev_time_min");

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         d_density_grad_tol = 
            db->getDoubleArray("d_density_grad_tol");
         d_density_grad_time_max = 
            db->getDoubleArray("d_density_grad_time_max");
         d_density_grad_time_min = 
            db->getDoubleArray("d_density_grad_time_min");

      }  else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         d_pressure_dev_tol = 
            db->getDoubleArray("d_pressure_dev_tol");
         d_pressure_dev = 
            db->getDoubleArray("d_pressure_dev");
         d_pressure_dev_time_max = 
            db->getDoubleArray("d_pressure_dev_time_max");
         d_pressure_dev_time_min = 
            db->getDoubleArray("d_pressure_dev_time_min");

      } else if  (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         d_pressure_grad_tol = 
            db->getDoubleArray("d_pressure_grad_tol");
         d_pressure_grad_time_max = 
            db->getDoubleArray("d_pressure_grad_time_max");
         d_pressure_grad_time_min = 
            db->getDoubleArray("d_pressure_grad_time_min");

      } 

   }

}

/*
*************************************************************************
* 从输入数据库读入Dirichlet边界条件.                                    *
*************************************************************************
*/

void Euler::readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                           string& db_name,
                                           int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
   assert(!db_name.empty());
#endif
#if (NDIM == 2)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_edge_density,
                      d_bdry_edge_velocity,
                      d_bdry_edge_pressure);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_density,
                      d_bdry_face_velocity,
                      d_bdry_face_pressure);
#endif
}

void Euler::readStateDataEntry(tbox::Pointer<tbox::Database> db,
                               const string& db_name,
                               int array_indx,
                               tbox::Array<double>& density,
                               tbox::Array<double>& velocity,
                               tbox::Array<double>& pressure)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
   assert(!db_name.empty());
   assert(array_indx >= 0);
   assert(density.getSize() > array_indx);
   assert(velocity.getSize() > array_indx*NDIM);
   assert(pressure.getSize() > array_indx);
#endif

   if (db->keyExists("density")) {
      density[array_indx] = db->getDouble("density");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`density' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("velocity")) {
      tbox::Array<double> tmp_vel(0);
      tmp_vel = db->getDoubleArray("velocity");
      if (tmp_vel.getSize() < NDIM) {
         TBOX_ERROR(d_object_name << ": "
                    << "Insufficient number `velocity' values"
                    << " given in " << db_name
                    << " input database." << endl);
      }
      for (int iv = 0; iv < NDIM; iv++) {
         velocity[array_indx*NDIM+iv] = tmp_vel[iv];
      }
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`velocity' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("pressure")) {
      pressure[array_indx] = db->getDouble("pressure");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`pressure' entry missing from " << db_name 
         << " input database. " << endl);
   }

}





