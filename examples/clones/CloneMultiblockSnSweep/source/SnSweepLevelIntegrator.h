 //
// 文件名:  SnSweepLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  Sn输运计算的网格层时间积分算法实现
//

#ifndef included_SnSweepLevelIntegrator
#define included_SnSweepLevelIntegrator
 
#include "JASMIN_config.h"

#include "SnSweep.h"
#include "TimeIntegratorLevelStrategy.h"

#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "ReductionIntegratorComponent.h"
#include "SweepingIntegratorComponent.h"
#include "CloneBcastComponent.h"
#include "CloneReductionComponent.h"

#define TESTING

using namespace JASMIN;

/**
 * @brief 该类以Sn角方向离散+隐式迎风差分格式+源迭代法,
 * 实现多群中子输运方程的网格层时间积分算法.
 *
 * 该类实现 advanceLevel() 函数时,
 * -# 利用数值构件实现外源项和散射源项的计算;
 * -# 利用扫描构件实现Sn扫描计算.
 *
 * 该类对象从输入文件读取如下参数:
 *    - dt \n
 *       双精度浮点型, 表示时间步长.
 *    - error_tole \n
 *       双精度浮点型, 源迭代的相对误差阈值.
 *    - max_sweep_iter \n 
 *       双精度浮点型, 源迭代的最大迭代次数.
 */
class SnSweepLevelIntegrator :
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    */     
   SnSweepLevelIntegrator(const string& object_name,
                   tbox::Pointer<tbox::Database> input_db,
                   SnSweep * snsweep_model,
                   tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geometry);

   /**
    * @brief 析构函数.
    */
   virtual ~SnSweepLevelIntegrator();

   /// @name 重载类 algs::TimeIntegratorLevelStrategy 的虚函数
   //@{

   /**
    * @brief 初始化网格层时间积分算法. 
    */
   void 
   initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 初始化指定网格层的数据片.
    *
    * 初值构件类 algs::InitializeIntegratorComponent 支撑该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向待初始化网格层所在的网格片层次结构.
    * @param level_number 输入参数, 整型, 表示待初始化的网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示该网格层可被进一步细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param old_level 输入参数, 指针, 指向初始化网格层的旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示初始化为数据片调度内存空间. 
    *
    */
   void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                       const int    level_number,
                       const double init_data_time,
                       const bool   can_be_refined,
                       const bool   initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level =
                               tbox::Pointer< hier::BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /**
    * @brief 获取指定网格层的时间步长. 
    *
    * 步长构件类 algs::DtIntegratorComponent
    * 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
    *
    * @return 双精度浮点型, 表示网格层的时间步长.
    */
   double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层积分一个时间步. 
    *
    * 数值构件类 algs::NumericalIntegratorComponent 
    * 和扫描构件类 algs::SweepingIntegratorCompnent
    * 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向待积分的网格层.
    * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
    * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
    * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
    * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
    * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
    * @param hierarchy_step_number, 输入参数, 整型, 表示网格片层次结构的积分步数,
    *                               也就是最粗网格层的积分步数.
    * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
    *
    * @return 整型, 表示该时间步积分的状态, 其值的含义由用户负责解释.
    */ 
   int
   advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time, 
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step, 
                const int    hierarchy_step_number,
                double&      actual_dt);

   /**
    * @brief 在指定的网格层上, 时间步积分完毕, 存储最新的获得数值解, 
    * 更新网格层的状态到新的时刻, 完成下一个时间步积分的准备工作.
    *
    * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
    * 
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    */
   void acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                  const double new_time, 
                  const bool deallocate_data);
   // @}

private:
   void getFromInput(tbox::Pointer<tbox::Database> db );

   // 对象名称.
   string d_object_name;

   // 网格片时间积分算法对象.
   SnSweep *d_snsweep_model;
   
   // 网格几何
   tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > d_grid_geometry;

   double d_dt;  // 时间步长
   double d_error_tole; // 误差阈值
   int    d_max_sweep_iter; // 最大迭代次数
   int    d_max_num_cells_swept;


   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >  d_init_intc;

   // 数值构件: 求外源项、散射源项和流通量。
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_sourcew_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_scatter_intc;
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_vflux_intc;
  
   // 扫描构件: 隐式迎风格式求解角通量.
   tbox::Pointer< algs::SweepingIntegratorComponent<NDIM> > d_solve_intc;

   // 规约构件: 计算源迭代误差.
   tbox::Pointer< algs::ReductionIntegratorComponent<NDIM> > d_iter_error_intc;

   // 复制构件: 接收数值解.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> > d_new_2_current_intc;

   // 内存构件 
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_new_context_intc;
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_scratch_context_intc;

   tbox::Pointer< algs::CloneBcastComponent<NDIM> > d_bcast_solution;
   tbox::Pointer< algs::CloneReductionComponent<NDIM> > d_collection_solution;

   tbox::Pointer<tbox::Timer> t_compute_sweep;
   tbox::Pointer<tbox::Timer> t_comm_bcast;

#ifdef TESTING
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_diff_intc;
#endif

};

#endif
