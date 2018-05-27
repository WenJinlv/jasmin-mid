 //
// 文件名:  MBLinAdvFederationLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  线性对流问题的联邦网格层时间积分算法的实现.
//


#include "MBLinAdvLevelIntegrator.h" 
#include "MBLinAdvFederationLevelIntegrator.h" 
#include "FederationComponentManager.h" 
#include "FederationMultiblockPatchHierarchy.h" 
#include "FederationMultiblockPatchLevel.h" 
#include "tbox/IEEE.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif


/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
MBLinAdvFederationLevelIntegrator::MBLinAdvFederationLevelIntegrator(
        const string& object_name,
        const int nfederals,
        MBLinAdvLevelIntegrator * level_integrator[],
        MBLinAdvFederation* patch_strategy )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(nfederals <= MAX_FEDERALS);
   for(int fd=0; fd<nfederals; fd++) {
       assert(level_integrator[fd]!=NULL);
   }
   assert(patch_strategy!=NULL);
#endif

   d_object_name = object_name;
   d_nfederals = nfederals;
   for(int fd=0; fd<d_nfederals; fd++) {
       d_level_integrator[fd] = level_integrator[fd];
   }
   d_patch_strategy = patch_strategy;
   
}

/*************************************************************************
 *
 * 析构函数.
 *
 ************************************************************************/
MBLinAdvFederationLevelIntegrator::~MBLinAdvFederationLevelIntegrator() {

}

/*************************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 * 该函数创建了1个内存构件，1个联邦构件.
 *************************************************************************/
void MBLinAdvFederationLevelIntegrator::initializeLevelIntegrator(
            tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

    // 联邦构件管理器.
    tbox::Pointer<algs::FederationComponentManager<NDIM> > 
                                        federation_manager = manager;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!federation_manager.isNull());
#endif

    // 初始化各个邦员.
    for(int fd=0; fd<d_nfederals; fd++) {
        d_level_integrator[fd]->initializeLevelIntegrator(
                                federation_manager->getFederalManager(fd));
    } 

    // 创建联邦构件.
    tbox::Pointer< algs::IntegratorComponentManager<NDIM> > fd_manager = 
                         federation_manager->getFederationManager();

    d_federation_intc = new algs::FederationIntegratorComponent<NDIM>(
                                           "FEDERATION",
                                           d_patch_strategy,
                                           fd_manager);

}


/*************************************************************************
 * 初始化指定网格层的数据片.
 ************************************************************************/
void MBLinAdvFederationLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     tbox::Pointer<hier::FederationMultiblockPatchHierarchy<NDIM> > fd_hier = hierarchy;
     tbox::Pointer<hier::FederationMultiblockPatchLevel<NDIM> > fd_old_level = old_level;
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!fd_hier.isNull());
     assert(d_nfederals==fd_hier->getNumberFederals());
     assert(d_nfederals==fd_hier->getNumberFederals());
     if(!fd_old_level.isNull()) 
        assert(d_nfederals==fd_old_level->getNumberFederals());
#endif

     for(int fd=0; fd<d_nfederals; fd++) {
         tbox::Pointer< hier::MultiblockPatchLevel<NDIM> > mb_level;
         if(!fd_old_level.isNull()) mb_level = 
                    fd_old_level->getMultiblockPatchLevelForFederal(fd);
         d_level_integrator[fd]->initializeLevelData(fd_hier->getHierarchy(fd),
                                    level_number,
                                    init_data_time,
                                    can_be_refined,
                                    initial_time,
                                    mb_level,
                                    allocate_data);
     }

}

/*************************************************************************
 * 计算时间步长.
 ************************************************************************/
double MBLinAdvFederationLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{

     tbox::Pointer<hier::FederationMultiblockPatchLevel<NDIM> > fd_level = level;
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!fd_level.isNull());
#endif

     double federal_dt=tbox::IEEE::getDBL_MAX();
     for(int fd=0; fd<d_nfederals; fd++) {
         double fdt = d_level_integrator[fd]->getLevelDt(
                              fd_level->getMultiblockPatchLevelForFederal(fd),
                              dt_time,
                              initial_time, 
                              flag_last_dt,
                              last_dt);
         if(fdt<federal_dt) federal_dt=fdt;
     }

     return(federal_dt);

}

/*************************************************************************
 * 向前积分一个时间步.
 ************************************************************************/
int MBLinAdvFederationLevelIntegrator::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt)
{
     tbox::Pointer<hier::FederationMultiblockPatchLevel<NDIM> > fd_level = level;
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!fd_level.isNull());
#endif
 
     // 开辟演算数据片的内存. 
     for(int fd=0; fd<d_nfederals; fd++) {
         d_level_integrator[fd]->allocateFederalXferData(
                 fd_level->getMultiblockPatchLevelForFederal(fd),
                 current_time);
     }

     // 交换联邦边界.
     d_federation_intc->transferLevelData(level, current_time);

     // 计算.
     for(int fd=0; fd<d_nfederals; fd++) {
         d_level_integrator[fd]->advanceLevel(
                                    fd_level->getMultiblockPatchLevelForFederal(fd),
                                    current_time,
                                    predict_dt,
                                    first_step);
     }
   
     // 释放演算内存空间.
     for(int fd=0; fd<d_nfederals; fd++) {
         d_level_integrator[fd]->deallocateFederalXferData(
                 fd_level->getMultiblockPatchLevelForFederal(fd));
     }

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
void MBLinAdvFederationLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{

     tbox::Pointer<hier::FederationMultiblockPatchLevel<NDIM> > fd_level = level;
#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!fd_level.isNull());
#endif

     for(int fd=0; fd<d_nfederals; fd++) {
         d_level_integrator[fd]->acceptTimeDependentSolution(
                                 fd_level->getMultiblockPatchLevelForFederal(fd),
                                 new_time,
                                 deallocate_data);
     }

}

