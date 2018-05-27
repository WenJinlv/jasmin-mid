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
        LinAdv* patch_strategy,
        const bool use_time_refinement)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(patch_strategy!=NULL);
#endif

   d_object_name = object_name;
   d_patch_strategy = patch_strategy;
   d_use_time_refinement = use_time_refinement;
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

    //数值构件: 标记待细化的网格单元.
    d_tag_intc  = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "TAG_CELLS",
                                                        d_patch_strategy,
                                                        manager);

    //同步构件: 沿粗细网格层交接面校正流. 
    d_sync_flux_intc = new algs::SynchronizeIntegratorComponent<NDIM>(
                                                       "SYNC_INTERFACE",
                                                       "INTERFACE_SYNC",
                                                        d_patch_strategy,
                                                        manager);

    //同步构件: 细网格层校正粗网格层.
    d_sync_overlay_intc = new algs::SynchronizeIntegratorComponent<NDIM>(
                                                       "SYNC_OVERLAY",
                                                       "OVERLAY_SYNC",
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

     // 存储粗细交界面的通量.
     if(level->getLevelNumber()>0) {
         d_sync_flux_intc -> storeInterfaceFluids(level,first_step);
     } 

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

/*************************************************************************
 *
 * 细网格层同步相邻的粗网格层.
 *
 * 该函数用细网格层数据矫正粗网格层的守恒量 u. 
 ************************************************************************/
void LinAdvLevelIntegrator::synchronizeCoarserLevel(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
      const int finer_level_number, 
      const double sync_time,
      const double coarser_old_time)
{ 

   tbox::Pointer< hier::PatchLevel<NDIM> > finer_level =
                 hierarchy -> getPatchLevel(finer_level_number);
   tbox::Pointer< hier::PatchLevel<NDIM> > coarser_level =
                 hierarchy -> getPatchLevel(finer_level_number-1);

   // 校正相邻粗网格层的通量. 
   d_sync_flux_intc->synchronizeCoarserLevel(finer_level,sync_time);

   // 粗网格层上重新计算守恒量 u.
   double dt = sync_time-coarser_old_time;
   d_diff_intc->computing(coarser_level,
                          coarser_old_time,
                          dt,
                          false);

   // 校正被细网格层覆盖的粗网格层计算区域.
   d_sync_overlay_intc->synchronizeCoarserLevel(finer_level, sync_time);

}

/*************************************************************************
 *
 * 标记待细化的网格单元。
 *
 ************************************************************************/
void LinAdvLevelIntegrator::tagCellsForRefinement(
               const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
               const int level_number,
               const double error_data_time,
               const int tag_index,
               const bool initial_time)
{
   tbox::Pointer< hier::BasePatchLevel<NDIM> > level =
                 hierarchy -> getPatchLevel(level_number);

   // 开辟数据片<uval,scratch>内存（标记细化单元时需要该数据片）
   d_alloc_scratch_data->allocatePatchData(level, error_data_time);

   // 设置存储标记值的数据片索引号.
   d_patch_strategy->setTagIndex(tag_index);

   // 标记细化单元.
   d_tag_intc->computing(level,
                         error_data_time,
                         0.0, true, initial_time);

   d_patch_strategy->clearTagIndex();

   // 释放内存.
   d_alloc_scratch_data->deallocatePatchData(level);

}

/*************************************************************************
 *
 * 返回时间积分的步进模式.
 *
 ************************************************************************/
bool LinAdvLevelIntegrator::usingRefinedTimestepping()
{
   return(d_use_time_refinement);
}


