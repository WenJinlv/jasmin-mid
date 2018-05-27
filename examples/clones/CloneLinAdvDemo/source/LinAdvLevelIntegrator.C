 //
// 文件名:  LinAdvLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  线性对流问题的网格层时间积分算法的实现.
//


#include "LinAdvLevelIntegrator.h" 
#include "CommunicatorDatabase.h" 

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
 * 该函数创建了内存构件，初值构件, 数值构件, 时间步长构件和复制构件各1个. 
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

    //初值构件: 管理当前值数据片的内存以及初始化
    d_init_set_value = new algs::InitializeIntegratorComponent<NDIM>(
                                                       "INIT_SET_VALUE",
                                                        d_patch_strategy,
                                                        false, manager);

    //步长构件: 计算稳定性时间步长. 
    d_step_size = new algs::DtIntegratorComponent<NDIM>("STEP_SIZE",
                                                   d_patch_strategy,
                                                   manager);

    //数值构件: 计算通量f 和守恒量u.
    d_core_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "COMPUTE",
                                                        d_patch_strategy,
                                                        manager);

    //克隆广播构件:将数值解从隆主网格层广播到克隆网格层.
    d_bcast_solution = new algs::CloneBcastComponent<NDIM>(
                                                       "BCAST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

    //克隆汇集构件:将数值解从克隆网格层汇集到属主网格层.
    d_reduction_solution = new algs::CloneReductionComponent<NDIM>(
                                                       "REDUCTION_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

    //克隆脉动构件:将数值解在克隆网格层之间脉动.
    d_pulsation_solution = new algs::ClonePulsationComponent<NDIM>(
                                                       "PULSATION_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

    //克隆汇集构件:初始化完毕，将数值解从克隆网格层汇集到属主网格层.
    d_initred_solution = new algs::CloneReductionComponent<NDIM>(
                                                       "INIT_COLL_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

    //数值构件: 初始化时刻, 在克隆网格层上，填充数据片影像区.
    d_init_fill_ghost = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "INIT_FILL_GHOST",
                                                        d_patch_strategy,
                                                        manager,
                                                        true);

    //数值构件: 初始化时刻, 后处理属主网格层的数据.
    d_initpost_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "INIT_POST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager,
                                                        true);

    //数值构件: 后处理属主网格层的数据.
    d_postprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "POST_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);


    //数值构件: 前处理克隆网格层的数据.
    d_prepprocess_computation = new algs::NumericalIntegratorComponent<NDIM>(
                                                       "PREP_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

    //复制构件: 属主网格层接受数值解.
    d_copy_solution = new algs::CopyIntegratorComponent<NDIM>(
                                                       "COPY_SOLUTION",
                                                        d_patch_strategy,
                                                        manager);

}


/*************************************************************************
 * 初始化指定网格层的数据片.
 *
 * 注解: 该函数调用了初值构件的成员函数initializeLevelData(),
 * 后者又进一步自动调用 d_patch_strategy->initializePatchData(), 
 * 完成数据片<uval,current>的初始化.
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
     // 在克隆网格层上计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {
        // 赋初值.
        d_init_set_value->initializeLevelData(hierarchy,
                                              level_number,
                                              init_data_time,
                                              initial_time,
                                              old_level);

       // 填充影像区.
       d_init_fill_ghost->computing(hierarchy->getPatchLevel(level_number), 
                                    init_data_time, 0, true, initial_time);

       // 汇集数值解到属主网格层.
       d_initred_solution->init_reduce(hierarchy->getPatchLevel(level_number),
                                         xfer::RefineUtilities::SUM);
    }
    hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

    // 在属主网格层上后处理数值解.
    d_initpost_computation->computing(hierarchy->getPatchLevel(level_number),
                                      init_data_time,
                                      0,
                                      true,
                                      initial_time);
}

/*************************************************************************
 * 计算时间步长.
 *
 * 注解: 该函数调用了步长构件的成员函数getLevelDt(),
 * 后者又进一步调用 d_patch_strategy->getPatchDt(), 
 * 逐个网格片地计算时间步长, 最后返回最小步长.
 ************************************************************************/
double LinAdvLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
     double step_size;   

     // 在克隆网格层上计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {
        step_size = d_step_size -> getLevelDt(level,
                                     dt_time,
                                     initial_time, 
                                     flag_last_dt,
                                     last_dt,
                                     false);
     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

     return(step_size);
}

/*************************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数首先调用内存构件为新值数据片申请内存空间；
 * 然后通过数值构件的成员函数的computing(),
 * 调用 d_patch_strategy->computeOnPatch(), 
 * 完成新时刻的通量 f 和守恒量 u 的计算.
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

     // 在克隆网格层上计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {

        // 申请内存
        if(first_step) { 
            d_alloc_new_data->allocatePatchData(level, current_time+predict_dt);
        }

        // 广播数值解.
        d_bcast_solution->bcast(level);

        // 前处理: 元素值*（clone_number+1)
        d_prepprocess_computation->computing(level, current_time, predict_dt);

        // 脉动计算: 元素值*（clone_number+1), rollback=true
        d_pulsation_solution->pulse(level,current_time,predict_dt,true);

        // 计算通量
        d_core_computation -> computing(level, current_time, predict_dt);

        // 汇集数值解.
        d_reduction_solution->reduce(level, xfer::RefineUtilities::SUM);

        actual_dt = predict_dt;

     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

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
     // 存储数值解.
     d_copy_solution ->copyPatchData(level, new_time);

     // 后处理数值解.
     d_postprocess_computation->computing(level,new_time,0);

     // 在克隆网格层上计算.
     hier::CommunicatorDatabase<NDIM>::getDatabase()->beginCloneComputing();
     {
        if(deallocate_data) d_alloc_new_data->deallocatePatchData(level);
     }
     hier::CommunicatorDatabase<NDIM>::getDatabase()->finalizeCloneComputing();

}


