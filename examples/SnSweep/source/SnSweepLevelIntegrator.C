 //
// 文件名:  SnSweepLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  Sn输运计算的网格层时间积分算法实现
//

#include "SnSweepLevelIntegrator.h" 
#include "tbox/IEEE.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif
#include <iomanip>
using namespace std;
using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
SnSweepLevelIntegrator::SnSweepLevelIntegrator(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        SnSweep *snsweep_model,
        tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geometry)
:
   d_object_name(object_name),
   d_snsweep_model(snsweep_model),
   d_grid_geometry(grid_geometry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(snsweep_model!=NULL);
#endif
   getFromInput(input_db);

#ifdef TESTING
   snsweep_model->setSolutionCoeff(d_dt);
#endif
}

/****************************************************************
 * 析构函数.
 ****************************************************************/
SnSweepLevelIntegrator::~SnSweepLevelIntegrator() {
}

/****************************************************************
 * 创建构件
 ****************************************************************/
void SnSweepLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{
  // 初值构件 : 初始化新创建的网格层.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_snsweep_model,
                                                              false, manager);

  // 数值构件: 求外源项、散射源项和流通量。
  d_sourcew_intc  = new algs::NumericalIntegratorComponent<NDIM>("SOURCEW",
                                                              d_snsweep_model,
                                                              manager);
  d_scatter_intc  = new algs::NumericalIntegratorComponent<NDIM>("SCATTER",
                                                              d_snsweep_model,
                                                              manager);
  d_vflux_intc    = new algs::NumericalIntegratorComponent<NDIM>("VFLUX",
                                                              d_snsweep_model,
                                                              manager);
#ifdef TESTING
  // 数值构件: 比对数值解和精确解。
  d_diff_intc    = new algs::NumericalIntegratorComponent<NDIM>("DIFF_RESULT",
                                                              d_snsweep_model,
                                                              manager);
#endif


  // 扫描构件: 隐式迎风格式求解角通量.
  d_solve_intc  = new algs::SweepingIntegratorComponent<NDIM>("SOLVE",
                                                              d_snsweep_model,
                                                              manager);
  d_solve_intc->setPerformanceParameters(d_max_num_cells_swept,"FILO",true);


  // 规约构件: 计算源迭代误差.
  d_iter_error_intc = new algs::ReductionIntegratorComponent<NDIM>("ITER_ERROR",
                                                              MPI_MAX,
                                                              d_snsweep_model,
                                                              manager);


  // 复制构件: 接收数值解.
  d_new_2_current_intc = new algs::CopyIntegratorComponent<NDIM>("NEW_2_CUR",
                                                              d_snsweep_model,
                                                              manager);

  // 内存构件: 当前值数据片, 新值数据片, 演算数据片.
  d_new_context_intc = new algs::MemoryIntegratorComponent<NDIM>("NEW",
                                                              d_snsweep_model,
                                                              manager);
  d_scratch_context_intc = new algs::MemoryIntegratorComponent<NDIM>("SCRATCH",
                                                              d_snsweep_model,
                                                              manager);

}


/****************************************************************
 * 初始化网格层数据.
 * 
 * 该函数按以下步骤执行:
 * (1) 为所有数据片调度内存空间;
 * (2) 初始化
 * (3) 为扫描构件计算Sn扫描过程中的有向数据依赖关系
 ****************************************************************/
void SnSweepLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data) 
{
     tbox::Pointer<hier::BasePatchLevel<NDIM> > level
           = hierarchy->getPatchLevel(level_number);

     // 按初始条件初始化角通量(fi,current)数据片.
     d_init_intc->initializeLevelData(hierarchy,
                                      level_number,
                                      init_data_time,
                                      initial_time,
                                      old_level,
                                      false);

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double SnSweepLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
     return(d_dt);   //采用固定的时间步长.
}

/****************************************************************
* 积分一个时间步长.
*
* 该函数按以下步骤执行:
* (1) 调用数值构件, 计算当前时刻的外源项.
* (2) 源迭代循环, 直到流通量的相对误差小于阈值, 或迭代次数过多:
*     (2.1) 调用数值构件, 计算散射源
*     (2.2) 调用扫描构件, 更新中子角通量(基于迎风格式).
*     (2.3) 计算流通量的相对迭代误差.
****************************************************************/
int SnSweepLevelIntegrator::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(predict_dt>=min_dt && predict_dt<=max_dt);
#endif

    NULL_USE(hierarchy_step_number);

    if(first_step) {
        // 为数据片调度内存空间.
        d_new_context_intc->allocatePatchData(level,current_time);
        d_scratch_context_intc->allocatePatchData(level,current_time);

        // 计算Sn扫描的有向数据依赖关系.
        d_solve_intc->setDataDependency(level,current_time,0.0);
    }

    actual_dt = predict_dt;
#ifdef TESTING
    d_snsweep_model->setTimeStepNumber(hierarchy_step_number);
#endif

    // 计算外源项
    d_sourcew_intc->computing(level,current_time,actual_dt);
    int iter;
    double wf_err_norm;
    for ( iter=0; iter<d_max_sweep_iter; iter++ ) {
         // 更新散射源
         d_scatter_intc->computing(level,current_time,actual_dt);

         // 更新角通量(基于隐式迎风格式).
         d_solve_intc->sweepingOnLevel(level,current_time,actual_dt); 

         // 更新流通量.
         d_vflux_intc->computing(level,current_time,actual_dt);

         // 计算流通量误差.
         wf_err_norm = 0.0;
         d_iter_error_intc->reduction(&wf_err_norm,1,level,current_time,actual_dt,false);

         if ( wf_err_norm < d_error_tole ) break;
    }

    tbox::pout<<std::setprecision(20)<<"#iter ="<<iter<<", error="<<wf_err_norm;

    // 设置新值数据片的时刻.
    d_new_context_intc->setTime(level,current_time+actual_dt);

#ifdef TESTING
    d_diff_intc->computing(level,current_time,actual_dt);
    tbox::pout<<"   Passed!";
#endif
    tbox::pout<<endl;

    return(0);

}

/****************************************************************
 * 接收数值解.
 *****************************************************************/
void SnSweepLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{
    d_new_2_current_intc->copyPatchData(level,new_time);
 
}

/********************************************************************
 *  从输入文件读取数据成员. 
 ********************************************************************/
void SnSweepLevelIntegrator::getFromInput(tbox::Pointer<tbox::Database> db )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif
   d_dt = db->getDouble("dt");
   d_error_tole = db->getDouble("error_tole");
   d_max_sweep_iter = db->getInteger("max_sweep_iter");
   d_max_num_cells_swept = db->getIntegerWithDefault(
                                     "max_num_cells_swept",
                                      tbox::IEEE::getINT_MAX());

}

