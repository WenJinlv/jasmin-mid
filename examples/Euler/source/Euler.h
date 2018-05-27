//
// 文件名:	Euler.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 170 $
// 修改  :	$Date: 2007-06-27 08:44:22 $
// 描述  :	Euler方程组：单网格片数值计算子程序.
//
 
#ifndef included_Euler
#define included_Euler

#ifndef included_JASMIN_config
#include "JASMIN_config.h"
#endif

#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Serializable.h"

#include "CellVariable.h"
#include "CellData.h"
#include "FaceVariable.h"
#include "FaceData.h"
#include "VariableContext.h"
#include "Patch.h"
#include "UniRectangularGridGeometry.h"
#include "RefinePatchStrategy.h"
#include "BoundaryUtilityStrategy.h"
#include "BoundaryDefines.h"
#include "StandardComponentPatchStrategy.h"
#include "IntegratorComponent.h"
#include "EulerEquationGodunov.h"

#include "JaVisDerivedDataStrategy.h"
#include "JaVisDataWriter.h"

using namespace JASMIN;


/**
 * @brief 类 Euler 提供单网格片上数值计算子程序, 求解Euler方程组.
 * 
 * 实现网格片时间积分算法基类 algs::StandardComponentPatchStrategy.
 * 在均匀矩形的单网格片上, 求解Euler方程组.
 *
 * 时间离散采用显式格式, 空间离散采用Goundov格式, 
 * 直接调用类appu::EulerEquationGodunov<DIM>中的数值通量, 格式的精度由输入文件指定.
 */

class Euler : 
   public tbox::Serializable,
   public algs::StandardComponentPatchStrategy<NDIM>,
   public appu::BoundaryUtilityStrategy,
   public appu::JaVisDerivedDataStrategy<NDIM>
{
public:

   /*!
    * @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示该对象的名字.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    * @param grid_geom   输入参数, 指针, 指向均匀矩形网格几何.
    * @note
    * 该函数执行以下操作:
    *  - 创建变量;
    *  - 初始化数据;
    *  - 将当前对象注册到重启动管理器;
    *  - 如果需断点续算, 则从重启动文件中读入数据; 
    *  - 从输入数据库中读入数据, 该操作将覆盖从重启动文件中读入的数据.
    * 
    */

   Euler(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< geom::UniRectangularGridGeometry<NDIM> > grid_geom);  

    /**
     * @brief 析构函数.
     */
   ~Euler();


  ///@name 重载类 algs::StandardComponentPatchStrategy 的虚函数
  // @{

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


   /**
    * @brief 获取细化插值算子的格式宽度.
    * @return 整型向量, 细化插值算子的格式宽度.
    * @note
    *  该函数是对抽象基类 xfer::RefinePatchStrategy 中声明的虚函数的实现. 
    */

   hier::IntVector<NDIM> getRefineOpStencilWidth() const;

   /**
    * @brief 网格细化通信的后处理函数. 
    * @param fine   输入参数, 网格片类，表示位于细网格层的目的网格片.
    * @param coarse 输入参数, 网格片类，表示位于粗网格层的源网格片.
    * @param fine_box 输入参数, Box类， 表示在细网格层的目的网格片中, 将被填充的区域.
    * @param ratio  输入参数, 整型向量类，表示相对于粗网格层, 细网格层的网格细化率.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    * @note
    *  该函数是对抽象基类 xfer::RefinePatchStrategy 中声明的虚函数的实现. 
    */
   virtual void postprocessRefineOperator(hier::Patch<NDIM>& fine,
                                  const hier::Patch<NDIM>& coarse,
                                  const hier::Box<NDIM>& fine_box,
                                  const hier::IntVector<NDIM>& ratio,
                                  const string& intc_name);

   /**
    * @brief 获得粗化插值算子的格式宽度.
    * @return 整型向量, 粗化插值算子的格式宽度.
    * @note
    *  该函数是对抽象基类 xfer::CoarsenPatchStrategy 中声明的虚函数的实现. 
    */

   hier::IntVector<NDIM> getCoarsenOpStencilWidth() const;

   /**
    * @brief 网格粗化通信的后处理函数.
    * @param coarse     输入参数, 网格片类, 表示位于粗网格层的目的网格片.
    * @param fine 	输入参数, 网格片类, 表示位于细网格层的源网格片.
    * @param coarse_box 输入参数, Box类, 表示在粗网格层的目的网格片中, 将被填充的区域.
    * @param ratio 	输入参数, 整型向量, 表示相对于粗网格层, 细网格层的网格细化率.
    * @param intc_name  输入参数, 字符串, 表示积分构件的名称.
    * @note
    *  该函数是对抽象基类 xfer::CoarsenPatchStrategy 中声明的虚函数的实现.
    */
   virtual void postprocessCoarsenOperator(hier::Patch<NDIM>& coarse,
                                   const hier::Patch<NDIM>& fine,
                                   const hier::Box<NDIM>& coarse_box,
                                   const hier::IntVector<NDIM>& ratio,
                                   const string& intc_name);
 
   // @}

   /// @name 重载边界工具策略类中的虚函数
    // @{

   /*!
    * @brief 从输入数据库读取物理量在指定物理边界处的值.
    * 
    * @param db 输入参数, 指针, 指向输入数据库
    * @param db_name  输入参数, 字符串, 表示输入数据库名.
    * @param bdry_location_index  输入参数, 整型, 指定物理边界位置.
    * @note
    * - 该函数是对策略类 appu::BoundaryUtilityStrategy 中声明的虚函数的实现.
    * - 对三维问题, 该函数只读取面型边界的数据; 对二维问题, 只读取棱型边界的数据.
    * - 参数 bdry_location_index 表示使用Dirichlet边界条件的面或者边
    *   的位置索引, 用于指明哪个方向的物理边界为Dirichlet边界, 并且从输入文件中
    *   获取该边界的数据.
    */
   void readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                       string& db_name,
                                       int bdry_location_index);

   // @}

   ///@name 重载类 tbox::Serializable 的虚函数
   // @{

   /*!
    * @brief 将数据成员输出到重启动数据库.
    * @param db 输入参数, 指针, 指向重启动数据库.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   // @}

   ///@name 可视化相关函数
   // @{
   
   /*!
    * @brief 注册JaVis数据输出器到JASMIN框架, 同时, 注册可视化数据片.
    * @param viz_writer 输入参数, 指针, 表示 JaVis 数据输出器.
    */
   void registerJaVisDataWriter(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer);

   /**
    * @brief 由基本量计算导出量并将其打包到缓冲区.
    * @param buffer 	输入参数, 指针, 指向双精度浮点型缓冲区.
    * @param patch 	输入参数, 网格片类, 网格片对象.
    * @param region 	输入参数, Box类, 待打包数据的索引区域.
    * @param variable_name 输入参数, 字符串, 导出量的名称.
    * @param depth_id 输入参数, 整型, 分量索引号. 对于标量, 索引号应为0; 
                        对于矢量, 索引号应在0和 DIM-1 之间; 对于张量, 
                        索引号应在0和 DIM*(DIM -1)之间. 对于数组量, 
                        索引号应在0 和(数组量分量总数-1)之间.
    * @return  逻辑型, 真值表明导出量在该网格片上已被准确地定义.
    * @note 
    *  - 该函数实现了基类appu::JaVisDerivedDataStrategy中的虚函数.
    *    用来计算每个网格片上的注册在JaVis输出器中的导出量.
    *  - 该函数将参数variable_name指定的绘图数据输出到参数buffer
    *    指定的缓冲区中, 该缓冲区位于网格片patch的区域region中.
    */
   bool packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch<NDIM>& patch,
      const hier::Box<NDIM>& region,
      const string& variable_name,
      int depth_id);

  // @}

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

   // 创建所有变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
   void registerModelVariables();

   // 计算通量.
   void computeFluxOnPatch(hier::Patch<NDIM>& patch,
                           const double time,
                           const double dt);
   /*
    * 校正物理量.
    */
   void conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const double dt);

   /**
    * @brief 在积分构件中, 为数据片填充物理边界条件.
    *
    * @param patch     输入参数, 网格片类, 表示网格片.
    * @param fill_time 输入参数, 双精度浮点型, 表示当前时刻.
    * @param ghost_width_to_fill 输入参数, 整型向量, 表示物理边界影像区被填充的宽度.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    */
   virtual void setPhysicalBoundaryConditions_1(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);

   virtual void setPhysicalBoundaryConditions_2(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);

   /*!
    * @brief 标记待细化的网格单元.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param tag_time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param tag_index 输入参数, 整型, 表示存储标记值的数据片.
    * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
    */
   virtual void tagCellsOnPatch(hier::Patch<NDIM>& patch,
                                const double  tag_time,
                                const int     tag_index,
                                const bool    initial_time);
                                
   /*!
    * @brief 标记待细化的网格单元.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param tag_time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param tag_index 输入参数, 整型, 表示存储标记值的数据片.
    * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
    */
   virtual void CJ_tagCellsOnPatch(hier::Patch<NDIM>& patch,
                                const double  tag_time,
                                const int     tag_index,
                                const bool    initial_time);

   /*
    * 重设物理边界,由函数computeFluxOnPatch_1()调用.
    */
   void boundaryReset(hier::Patch<NDIM>& patch,
                      pdat::FaceData<NDIM,double>& traced_left,
                      pdat::FaceData<NDIM,double>& traced_right) const;
   
   /*
    * 读入边界条件数据,由函数getFromInput()调用.
    */
   void readStateDataEntry(tbox::Pointer<tbox::Database> db,
                           const string& db_name,
                           int array_indx,
                           tbox::Array<double>& density,
                           tbox::Array<double>& velocity,
                           tbox::Array<double>& pressure);

   /*
    * 两种三维校正通量的计算方法. 由函数updateFluxesOnPatch()调用.
    */
#if (NDIM==3)

   void compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                            const double dt);
   void compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                            const double dt);
#endif
   
   /*
    * 对象名.
    */
   string d_object_name;     

   /*
    * 网格几何、可视化数据输出器和绘图上下文 .
    */
   tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > d_grid_geometry;
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > d_javis_writer;
   tbox::Pointer<hier::VariableContext> d_plot_context;  

   tbox::Pointer<appu::EulerEquationGodunov<NDIM> > d_godunov_flux;

   /*
    * 存储数值解(密度, 速度和压力)及通量.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_density;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_velocity;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_pressure;
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_flux;
   /*
    * 燃烧函数.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_df;

   // 上下文: 当前值, 新值, 可视化.
   tbox::Pointer< hier::VariableContext > d_current;
   tbox::Pointer< hier::VariableContext > d_new;
   tbox::Pointer< hier::VariableContext > d_scratch;


   // 数据片索引号.
   int d_density_current_id, d_density_new_id, d_density_scratch_id, d_flux_new_id;
   int d_velocity_current_id, d_velocity_new_id, d_velocity_scratch_id;
   int d_pressure_current_id, d_pressure_new_id, d_pressure_scratch_id;
   int d_df_current_id, d_df_new_id, d_df_scratch_id;
   
   /*
    * 比热比.
    */
   double d_gamma;  
                  
   /*
    *  数值子程序中使用的参数.
    */
   string d_riemann_solve;
   int d_riemann_solve_int;               // 通量计算中使用的Riemann求解器
   int d_godunov_order;                   // Godunov格式的阶数
   string d_corner_transport;             // 三维有限差分中的通量校正方式
   hier::IntVector<NDIM> d_nghosts;       // 中心量变量的影像区宽度
   hier::IntVector<NDIM> d_fluxghosts;    // 通量(面心量)变量的影像区宽度

   string d_data_problem;
   int d_data_problem_int;                // 计算物理模型
   
   /*
    * 存储两点起爆模型问题的参数.
    */ 
   string d_model_problem;
   int    d_model_type;
   double d_left_point[NDIM];
   double d_right_point[NDIM];

   string d_detonation_model;
   double d_DJ_density;
   double d_CJ_velocity;
   double d_CJ_density;
   double d_CJ_pressure;
   double d_CJ_rb;      // CJ 参数rb, 见李德元老师书P97.
   double d_CJ_nb;      // CJ 参数nb, 见李德元老师书P97.

   double d_CJ_eng;  // 化学能, 初始内能.

   /*
    * SPHERE问题输入参数.
    */ 
   double d_radius;                       // 球半径
   double d_center[NDIM];                 // 球心坐标
   double d_density_inside;               // 球内密度
   double d_velocity_inside[NDIM];        // 球内速度
   double d_pressure_inside;              // 球内压力
   double d_density_outside;              // 球外密度
   double d_velocity_outside[NDIM];       // 球内速度
   double d_pressure_outside;             // 球内压力

   /*
    *  STEP问题的输入参数.
    */
   int d_number_of_intervals;               //计算区域个数
   tbox::Array<double> d_front_position;    //台阶位置
   tbox::Array<double> d_interval_density;  //台阶上方气体密度
   tbox::Array<double> d_interval_velocity; //台阶上方气体速度
   tbox::Array<double> d_interval_pressure; //台阶上方气体压力

   /*
    * 存储结点和棱(面)型边界的变量. 
    */
   tbox::Array<int> d_master_bdry_edge_conds;  //棱型边界
   tbox::Array<int> d_master_bdry_node_conds;  //结点型边界
#if (NDIM == 3)
   tbox::Array<int> d_master_bdry_face_conds;  //面型边界
#endif

   /*
    * 存储不同边界的标量和向量型的变量.
    */
   tbox::Array<int> d_scalar_bdry_edge_conds;  //棱型边界标量型变量
   tbox::Array<int> d_vector_bdry_edge_conds;  //棱型边界向量型变量
   tbox::Array<int> d_scalar_bdry_node_conds;  //结点型边界标量型变量
   tbox::Array<int> d_vector_bdry_node_conds;  //结点型边界向量型变量
#if (NDIM == 3)
   tbox::Array<int> d_scalar_bdry_face_conds;  //面型边界标量型变量
   tbox::Array<int> d_vector_bdry_face_conds;  //面型边界向量型变量
#endif

   /*
    * 存储不同类型相邻的边界的变量.
    */
#if (NDIM ==2)
   tbox::Array<int> d_node_bdry_edge;     //结点型边界相邻的棱型边界
#endif
#if (NDIM == 3)
   tbox::Array<int> d_edge_bdry_face;     //棱型边界相邻的面型边界
   tbox::Array<int> d_node_bdry_face;     //结点型边界相邻的面型边界
#endif


   /*
    * 存储面(三维)或棱(二维)上的DIRICHLET边界条件的变量.
    */
#if (NDIM ==2)
   tbox::Array<double> d_bdry_edge_density;
   tbox::Array<double> d_bdry_edge_velocity;
   tbox::Array<double> d_bdry_edge_pressure;
#endif
#if (NDIM ==3)
   tbox::Array<double> d_bdry_face_density;
   tbox::Array<double> d_bdry_face_velocity;
   tbox::Array<double> d_bdry_face_pressure;
#endif

   /*
    * 用于梯度检测法的网格细化标准.
    */
   tbox::Array<string> d_refinement_criteria;
   tbox::Array<double> d_density_dev_tol;
   tbox::Array<double> d_density_dev;
   tbox::Array<double> d_density_dev_time_max;
   tbox::Array<double> d_density_dev_time_min;
   tbox::Array<double> d_density_grad_tol;
   tbox::Array<double> d_density_grad_time_max;
   tbox::Array<double> d_density_grad_time_min;
   tbox::Array<double> d_density_shock_onset;
   tbox::Array<double> d_density_shock_tol;
   tbox::Array<double> d_density_shock_time_max;
   tbox::Array<double> d_density_shock_time_min;
   tbox::Array<double> d_pressure_dev_tol;
   tbox::Array<double> d_pressure_dev;
   tbox::Array<double> d_pressure_dev_time_max;
   tbox::Array<double> d_pressure_dev_time_min;
   tbox::Array<double> d_pressure_grad_tol;
   tbox::Array<double> d_pressure_grad_time_max;
   tbox::Array<double> d_pressure_grad_time_min;
   tbox::Array<double> d_pressure_shock_onset;
   tbox::Array<double> d_pressure_shock_tol;
   tbox::Array<double> d_pressure_shock_time_max;
   tbox::Array<double> d_pressure_shock_time_min;


};

#endif
