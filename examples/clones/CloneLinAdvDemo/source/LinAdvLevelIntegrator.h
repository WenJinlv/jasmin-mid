//
// 文件名:	LinAdvLevelIntegrator.h
// 软件包:	JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 250 $
// 修改  :	$Date: 2007-05-24 10:56:58$
// 描述  :	线性对流问题的网格层时间积分算法
//

#ifndef included_LinAdvLevelIntegrator
#define included_LinAdvLevelIntegrator
 
#include "TimeIntegratorLevelStrategy.h"
#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "CloneBcastComponent.h"
#include "CloneReductionComponent.h"
#include "ClonePulsationComponent.h"
using namespace JASMIN;

#include "LinAdv.h"

/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生,
 * 实现线性对流方程在网格层上的求解流程.
 */
class LinAdvLevelIntegrator 
: public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    * @param object_name     输入参数, 字符串, 表示对象名称.
    * @param patch_strategy  输入参数, 指针, 线性对流方程网格片积分算法.
    */     
   LinAdvLevelIntegrator(const string& object_name, LinAdv* patch_strategy);

   /**
    * @brief 析构函数.
    */
   virtual ~LinAdvLevelIntegrator();

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
    * 具体地，该函数完成以下操作：
    * - 若输入参数 initial_time 为真，则根据初始条件，
    *   为当前网格层上的所有数据片 <uval,current> 赋初值;
    * - 若输入参数 initial_time 为假（此时刚完成负载调整），
    *   则将数据片 <uval,current> 的值从旧网格层复制到新网格层.
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
    * 该函数调用了1个初值构件对象，该对象又进一步自动调用函数
    * LinAdv::initializePatchData(), 逐个网格片地完成初始时刻的数据初始化.
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
    * LinAdv::getPatchDt(), 逐个网格片地计算稳定时间步长.
    */
   double getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层向前积分一个时间步. 
    *
    * 具体地, 该函数完成以下操作：
    * -# 计算新时刻的通量@f$ f @f$, 存储到数据片<flux,new>中;
    * -# 计算新时刻的守恒量@f$ u @f$, 存储到数据片<uval,new>中;
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
    * @return 整型, 表示该时间步积分的状态. 
    *
    * @note
    * 该函数调用了2个数值构件对象，
    * 分别完成新时刻的通量 f 和守恒量 u 的计算.
    * 相关的网格片计算函数见 LinAdv::computeOnPatch().
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

   /*!@brief 线性对流的网格片积分算法. */
   LinAdv* d_patch_strategy;

   /*!@brief 内存构件: 管理新值数据片 */
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> > d_alloc_new_data;

   /*!@brief 初值构件: 管理u的当前值数据片的内存以及初始化 */
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> > d_init_set_value;

   /*!@brief 步长构件: 计算稳定性时间步长. */
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> > d_step_size;

   /*!@brief 数值构件: 计算通量 f 和物理量 u */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_core_computation;

   /*!@brief 数值构件: 在克隆网格层上，数值积分过程中前处理数值解. */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_prepprocess_computation;

   /*!@brief 数值构件: 在克隆网格层上，在初始时刻填充影像区. */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_init_fill_ghost;

   /*!@brief 克隆规约构件: 汇集数值解到属主网格层. */
   tbox::Pointer< algs::CloneReductionComponent<NDIM> > d_reduction_solution;

   /*!@brief 克隆脉动构件: 在克隆网格层之间脉动数据片. */
   tbox::Pointer< algs::ClonePulsationComponent<NDIM> > d_pulsation_solution;

   /*!@brief 数值构件: 初始化时刻，汇集数值解到属主网格层. */
   tbox::Pointer< algs::CloneReductionComponent<NDIM> > d_initred_solution;

   /*!@brief 克隆广播构件: 时间积分过程中，广播数值解到克隆网格层. */
   tbox::Pointer< algs::CloneBcastComponent<NDIM> > d_bcast_solution;

   /*!@brief 复制构件: 接受数值解. */
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> > d_copy_solution;

   /*!@brief 数值构件: 在属主网格层上，初始化时刻后处理数值解. */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_initpost_computation;

   /*!@brief 数值构件: 在属主网格层上，数值积分过程中后处理数值解. */
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> > d_postprocess_computation;

};

#endif
