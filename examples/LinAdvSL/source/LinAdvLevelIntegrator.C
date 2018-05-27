 //
// 文件名:  LinAdvLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  线性对流问题的网格层时间积分算法的实现.
//


#include "LinAdvLevelIntegrator.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif


/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
LinAdvLevelIntegrator::LinAdvLevelIntegrator(
        const string& object_name,
        LinAdv* patch_strategy )
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
LinAdvLevelIntegrator::~LinAdvLevelIntegrator() {

}

/*************************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 *
 * 该函数创建了2个内存构件，1个初值构件, 2个数值构件, 
 * 1个时间步长构件和1个复制构件. 
 * 这些构件所操作的数据片, 
 * 由函数 d_patch_strategy->initializeComponent() 指定.
 *
 *************************************************************************/
void LinAdvLevelIntegrator::initializeLevelIntegrator(
            tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

    //内存构件: 管理新值数据片 
    d_alloc_new_data = new algs::MemoryIntegratorComponent<NDIM>(
                                                       "ALLOC_NEW_DATA",
                                                        d_patch_strategy,
                                                        manager);

    // 内存构件: 管理演算数据片 
    d_alloc_scratch_data = new algs::MemoryIntegratorComponent<NDIM>(
                                                       "ALLOC_SCRATCH_DATA",
                                                        d_patch_strategy,
                                                        manager);

    //初值构件: 管理u的当前值数据片的内存以及初始化 
    d_init_set_value = new algs::InitializeIntegratorComponent<NDIM>(
                                                       "INIT_SET_VALUE",
                                                        d_patch_strategy,
                                                        false, manager);

    //步长构件: 计算稳定性时间步长. 
    d_step_size = new algs::DtIntegratorComponent<NDIM>("STEP_SIZE",
                                                        d_patch_strategy,
                                                        manager);

    //数值构件: 计算通量f 
    d_compute_flux = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "COMPUTE_FLUX",
                                                        d_patch_strategy,
                                                        manager);

    //数值构件: 更新守恒量u. 
    d_diff_intc  = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "CONSER_DIFF",
                                                        d_patch_strategy,
                                                        manager);

    //复制构件: 接受数值解. 
    d_copy_solution = new algs::CopyIntegratorComponent<NDIM>(
                                                       "COPY_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);
}


/*************************************************************************
 * 初始化指定网格层的数据片.
 *
 * 注解: 该函数调用了初值构件（d_init_set_value)，
 * 该构件又进一步自动调用 d_patch_strategy->initializePatchData(), 
 * 完成数据片<uval,current>的初始化.
 * 如果当前刚调整完负载，形成了新的网格层，那么该构件
 * 会将数据片<uval,current>的值从旧网格层复制到新网格层.
 *
 ************************************************************************/
void LinAdvLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     d_init_set_value->initializeLevelData(hierarchy,
                                           level_number,
                                           init_data_time,
                                           initial_time,
                                           old_level);
}

/*************************************************************************
 * 计算时间步长.
 *
 * 注解: 该函数调用了步长构件(d_step_size)，
 * 该构件对象又进一步调用 d_patch_strategy->getPatchDt(), 
 * 逐个网格片地计算时间步长.
 *
 ************************************************************************/
double LinAdvLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
 
     return(d_step_size -> getLevelDt(level,
                                      dt_time,
                                      initial_time, 
                                      flag_last_dt,
                                      last_dt,
                                      false));

}

/*************************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数调用了2个数值构件对象的computing()函数，
 * 分别完成新时刻的通量 f 和守恒量 u 的计算.
 *
 ************************************************************************/
int LinAdvLevelIntegrator::advanceLevel(
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
     assert(!level.isNull());
#endif

     // 开辟新值数据片的内存
     if(first_step) { 
         d_alloc_new_data ->allocatePatchData(level, current_time+predict_dt);
     }

     // 开辟演算数据片的内存
     d_alloc_scratch_data ->allocatePatchData(level, current_time+predict_dt);

     // 计算通量
     d_compute_flux -> computing(level,
                                  current_time,
                                  predict_dt);

     // 守恒差分.
     d_diff_intc->computing(level,
                            current_time,
                            predict_dt,
                            false);

     // 设置新值数据片的时戳
     d_alloc_new_data ->setTime(level, current_time+predict_dt);

     // 释放临时的内存空间.
     d_alloc_scratch_data ->deallocatePatchData(level);

     actual_dt = predict_dt;

     return(0); 
}

/*************************************************************************
 * 接收数值解.
 *
 * 注解: 该函数调用复制构件，
 * 将数据片<uval,new>的值复制到数据片<uval,current>中.
 *
 ************************************************************************/
void LinAdvLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{

   d_copy_solution ->copyPatchData(level, new_time);

   if(deallocate_data) d_alloc_new_data -> deallocatePatchData(level);

}


