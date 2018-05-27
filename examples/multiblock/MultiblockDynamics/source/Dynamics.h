//
// 文件  :	Dynamics.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 1.12 $
// 修改  :	$Date: 2006/08/07 03:13:21 $
// 描述  :	求解流体力学方程的用户网格片子程序.
//
 
#ifndef included_Dynamics
#define included_Dynamics

#include "StandardComponentPatchStrategy.h"
#include "MultiblockDeformingGridGeometry.h"
#include "DeformingGridInputUtilities2.h"
#include "CellVariable.h"
#include "NodeVariable.h"
#include "SideVariable.h"
#include "NodeData.h"
#include "JaVisDataWriter.h"
using namespace JASMIN;

/**
 * @brief 类 Dynamics 具体实现了求解流体力学方程的网格片子程序.
 *
 * 当前,该类基于IGA格式求解动量守恒方程.
 *
 * 该类的基本变量包括: 网格坐标和速度(结点量),压力、密度和能量(中心量).
 *
 * 该类需要从输入数据库读取如下参数. 
 * 
 *     - \b model \n
 *          字符串型, 模型名称. 
 *
 *     - \b a0 \n
 *          双精度浮点数, 一次黏性系数.
 *
 *     - \b b0 \n
 *          双精度浮点数, 二次黏性系数.
 *
 *     - \b initial_gamma \n
 *          双精度浮点数, 理想气体初始比热比.
 *
 *     - \b initial_density \n
 *          双精度浮点数, 初始密度.
 *
 *     - \b initial_energy \n
 *          双精度浮点数, 初始能量.
 *
 *     - \b initial_center_energy \n
 *          双精度浮点数, 初始中心能量.
 *
 *     - \b initial_velocity \n
 *          双精度浮点数, 初始速度.
 *
 *     - \b anc \n
 *          双精度浮点数, ANC方法系数.
 *
 *     - \b xi \n
 *          双精度浮点数, ANC方法系数.
 *  
 */
class Dynamics : 
   public algs::StandardComponentPatchStrategy<NDIM>,
   public tbox::Serializable
{
public:
   /**
    * 构造函数.
    */
   Dynamics(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer<appu::DeformingGridInputUtilities2> grid_tool,
         tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom);  

    /**
     * 析构函数.
     */
   ~Dynamics();
 
   /*!
    * @brief 初始化积分构件.
    *
    * 向指定的积分构件注册待操作的数据片.
    *
    * @param intc 输入参数, 指针, 指向待初始化的积分构件对象.
    */
    void initializeComponent( algs::IntegratorComponent<NDIM>* intc) const;

   /**
    * @brief 初始化数据片.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param intc_name 输入参数, 字符串, 表示初值构件的名称.
    */
    void initializePatchData(
      hier::Patch<NDIM>& patch,
      const double  time,
      const bool    initial_time,
      const string& intc_name);

   /*!
    * @brief 计算时间步长.
    *
    * 计算并返回指定网格片上的稳定时间步长.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt   输入参数, 整型, 表示前次积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示前次积分的时间步长.
    * @param intc_name 输入参数, 字符串, 表示时间步长构件的名称.
    *
    * @return 双精度浮点型, 表示时间步长.
    */
    double getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name);
   
   /*!
    * @brief 为指定构件完成数值计算.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示当前时刻.
    * @param dt    输入参数, 双精度浮点型, 表示时间步长.
    * @param initial_time  输入参数, 逻辑型, 表示当前是否为初始时刻.
    * @param intc_name 输入参数, 字符串, 表示数值构件的名称.
    */
    void computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name);

   /**
    * @brief 填充物理边界影像区.
    *
    * @param patch     输入参数, 网格片类, 表示网格片.
    * @param fill_time 输入参数, 双精度浮点型, 表示当前时刻.
    * @param ghost_width_to_fill 输入参数, 整型向量, 表示物理边界影像区被填充的宽度.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    */
    void setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);

   /**
    * 将Dynamics类对象的数据成员输出到指定的重启动数据库.
    *
    * 该函数是对类: tbox::Serializable 中定义的函数接口的实现.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);
 
   /**
    * 注册绘图量到JaVis数据输出器.           
    */
   void registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer);

   /**
    * 输出当前对象中所有数据成员.
    */
   void printClassData(ostream& os) const;

private:
   /**
    * 注册变量.
    */
   void registerModelVariables();

   /**
    * 预处理流体的基本量,包括压力和黏性等.
    */ 
   virtual void preprocessStateOnPatch(hier::Patch<NDIM>& patch, 
                             const double time, 
                             const double dt);

   /**
    * 计算网格边的压力.
    */ 
   virtual void computeSidePressureOnPatch(hier::Patch<NDIM>& patch, 
                             const double time, 
                             const double dt);

   /**
    * 计算网格结点加速度.
    */ 
    void computeSpeedupOnPatch(hier::Patch<NDIM>& patch);
    void computeSpeedupOnBox(hier::Patch<NDIM>& patch, const hier::Box<NDIM> &box); 

   /**
    * ANC方法调整结点速度.
    */ 
   void ANCOnPatch(hier::Patch<NDIM>& patch);
   void ANCOnBox(hier::Patch<NDIM>& patch,
	         const hier::Box<NDIM>& bSing,
		 tbox::Pointer<pdat::NodeData<NDIM,double> > velocity_tmp);

   /**
    * 计算网格结点速度.
    */ 
   void computeVelocityOnPatch(hier::Patch<NDIM>& patch, 
                             const double time, 
                             const double dt);
   /**
    * 后处理网格结点的速度:
    * 修正物理边界结点的速度.
    */ 
   void postprocessVelocityOnPatch(hier::Patch<NDIM>& patch, 
                             const double time, 
                             const double dt);

   /**
    * 计算网格结点的坐标、单元密度和能量.
    */
   void computeGridDensityEnergyOnPatch(hier::Patch<NDIM>& patch,
                            const double time,
                            const double dt);

   /*
    * 从输入数据库读入数据.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);
   /*
    * 从重启动数据库读入数据.
    */
   void getFromRestart();

   /*
    * 对象名.
    */
   string d_object_name;

   tbox::Pointer< geom::DeformingGridGeometry<NDIM> > d_grid_geometry;
   tbox::Pointer< appu::DeformingGridInputUtilities2 > d_grid_tool;

  /*
   * 动力学基本物理量,定义在二维变形网格的结点,包括网格坐标和流体速度.
   */
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_coordinates;   
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_velocity;

   /*
    * 热力学基本物理量,定义在二维变形网格的中心,包括流体密度和温度.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_density;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_energy;

   /*
    * 流体力学导出物理量,定义在二维变形网格的中心,包括流体压力和人为粘力.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_pressure;   // 总压力
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_side_pressure;   // 边压力
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_mass;       // 总质量
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_speedup;    // 加速度
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_viscosity;  // 人为粘性
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_paq;        // 压力和人为粘力之和
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_gamma;   

   /*
    * 描述数值解的辅助物理量,定义在二维变形网格的结点和中心.
    */
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_node_temp_var;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_cell_temp_var;

   tbox::Pointer< hier::VariableContext > d_new, d_current, d_scratch;

   int d_coordinates_current_id,
       d_coordinates_new_id,
       d_coordinates_scratch_id;

   int d_density_current_id,
       d_density_new_id,
       d_density_scratch_id;

   int d_velocity_current_id,
       d_velocity_new_id,
       d_velocity_scratch_id;

   int d_energy_current_id,
       d_energy_new_id;

   int d_pressure_scratch_id;
   int d_side_pressure_scratch_id;

   int d_viscosity_scratch_id,
       d_paq_scratch_id,
       d_mass_scratch_id,
       d_speedup_scratch_id;


   int d_gamma_current_id;

   /*
    * 以下为一些具体的物理参数和控制参数,需要从输入数据库或者重启动数据库中读入.
    */

    // 修正流体速度的控制物理量,省缺值为0.
    double d_anc;
    double d_xi;

    //二次和一次黏性系数
    double d_a0;
    double d_b0;

    /* 
     * 物理模型.
     */
   string d_model;

   /* 初始物理状态量 */
   double d_initial_velocity[NDIM]; // 初始流体速度
   double d_initial_density;        // 初始介质密度
   double d_initial_energy;        // 初始介质密度
   double d_initial_center_energy;        // 初始介质密度
   double d_initial_gamma;

};

#endif
