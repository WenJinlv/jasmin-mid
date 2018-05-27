 //
// 文件名:	FDMBSingDemoFederationLevelIntegrator.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 250 $
// 修改  :	$Date: 2007-05-24 10:56:58$
// 描述  :	联邦多块变形结构网格应用的网格层时间积分算法
//

#ifndef included_FDMBSingDemoFederationLevelIntegrator
#define included_FDMBSingDemoFederationLevelIntegrator
 
#include "FederationIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"

#include "FDMBSingDemoFederation.h"
#include "FDMBSingDemoLevelIntegrator.h"

using namespace JASMIN;

#define MAX_FEDERALS 3

/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生,
 * 实现网格层上的验证算法流程.
 */
class FDMBSingDemoFederationLevelIntegrator 
: public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    * @param object_name      输入参数, 字符串, 表示对象名称.
    * @param nfederals        输入参数, 整型, 表示邦员总数.
    * @param level_integrator 输入参数, 指针数组(长度为nfederals), 
    *                         表示各个邦员的网格层时间时间积分算法.
    * @param patch_strategy   输入参数, 指针, 线性对流方程联邦网格片积分算法.
    */     
   FDMBSingDemoFederationLevelIntegrator(
        const string& object_name,
        const int nfederals,
        FDMBSingDemoLevelIntegrator * level_integrator[],
        FDMBSingDemoFederation* patch_strategy );

   /**
    * @brief 析构函数.
    */
   virtual ~FDMBSingDemoFederationLevelIntegrator();

   ///@name 重载基类algs::TimeIntegratorLevelStrategy<NDIM>的函数
   //@{
   /**
    * @brief 初始化该积分算法: 创建所有计算需要的积分构件.
    * @param manager 输入参数, 指针, 指向积分构件管理器.
    */
   void initializeLevelIntegrator(
           tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 初始化指定网格层的数据片.
    *
    * @param hierarchy 输入参数, 指针, 指向待初始化网格层所在的网格片层次结构.
    * @param level_number 输入参数, 整型, 表示待初始化的网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示该网格层可被进一步细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param old_level 输入参数, 指针, 指向初始化网格层的旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示初始化为数据片调度内存空间. 
    *
    * @note
    * 该函数调用邦员的类似函数，初始化数据片。
    */
   void initializeLevelData(
                  const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                  const int    level_number,
                  const double init_data_time,
                  const bool   can_be_refined,
                  const bool   initial_time,
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
                  const bool allocate_data );

   /**
    * @brief 返回指定网格层的时间步长. 
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
    * @return 双精度浮点型, 表示网格层的时间步长.
    *
    * @note
    * 该函数调用了1个时间步长构件对象，该对象又进一步自动调用函数
    * FDMBSingDemo::getPatchDt(), 逐个网格片地计算稳定时间步长.
    */
   double getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层向前积分一个时间步. 
    * @param level 输入参数, 指针, 指向待积分的网格层.
    * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
    * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
    * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
    * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
    * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
    * @param hierarchy_step_number, 输入参数, 整型, 表示网格片层次结构的积分步数,
    *                               也就是最粗网格层的积分步数.
    * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
    * @return 整型, 表示该时间步积分的状态. 
    *
    * 该函数调用邦员的网格层时间积分函数，完成一个时间步的积分.
    */ 
   int advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time, 
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step, 
                const int    hierarchy_step_number,
                double&      actual_dt);

   /**
    * @brief 更新网格层的状态到新的时刻.
    * 
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    * 
    * @note
    * 该函数调用复制构件, 将数据从新值复制到当前值上下文的数据片.
    */
   void acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                  const double new_time, 
                  const bool deallocate_data);

   //@}
   
private:

   /*!@brief 对象名称. */
   string d_object_name; 

   /*!@brief 邦员总数. */
   int d_nfederals;

   /*!@brief 网格层时间积分算法. */
   FDMBSingDemoLevelIntegrator * d_level_integrator[MAX_FEDERALS];

   /*!@brief 网格片时间积分算法. */
   FDMBSingDemoFederation* d_patch_strategy;

   /*!@brief 联邦积分构件，内存构件 */
   tbox::Pointer<algs::FederationIntegratorComponent<NDIM> > d_federation_intc;

};

#endif
