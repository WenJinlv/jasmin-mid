//
// 文件名:  MDLevelIntegrator.C
// 软件包:  JASMIN applications
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 250 $
// 修改  :  $Date: 2007-05-24 10:56:58$
// 描述  :  网格层时间积分算法的实现.
//

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "MDLevelIntegrator.h" 
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ParticleData.h"


#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
MDLevelIntegrator::MDLevelIntegrator(
    const string& object_name,
    tbox::Pointer<tbox::Database> input_db,
    algs::StandardComponentPatchStrategy<NDIM>* patch_strategy,
    MD* md_model,
    tbox::Pointer< hier::GridGeometry<NDIM> > grid_geometry)
        :
        d_object_name(object_name),
        d_patch_strategy(patch_strategy),
        d_grid_geometry(grid_geometry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!grid_geometry.isNull());
    assert(patch_strategy != ((algs::StandardComponentPatchStrategy<NDIM>*)NULL));
#endif

    d_md_model = md_model ;

    bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if(from_restart) getFromRestart();
    getFromInput(input_db, from_restart);

}

/****************************************************************
* 析构函数.
*****************************************************************/
MDLevelIntegrator::~MDLevelIntegrator() {
}

/****************************************************************
* 创建并初始化所有积分构件.
*****************************************************************/
void MDLevelIntegrator::initializeLevelIntegrator(
            tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{

    // 初值构件 : 初始化新创建的网格层.
    d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                        d_patch_strategy,
                                                        false, manager);

    // 步长构件 : 计算时间步长.
    d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                        d_patch_strategy,
                                                        manager);

    // 粒子量通信构件: 负责粒子量的通信.
    d_particle_comm_intc = new algs::ParticleCommComponent<NDIM>("PARTICLE_COMM",
                                                        d_patch_strategy,
                                                        manager);

    // 数值构件: 更新粒子位置.
    d_particle_comp_intc = new algs::NumericalIntegratorComponent<NDIM>("PARTICLE_COMP",
                                                        d_patch_strategy,
                                                        manager);

    // 规约构件：计算粒子总数目.
    d_particle_reduce_intc = new algs::ReductionIntegratorComponent<NDIM>("PARTICLE_REDUCE",
                                                        MPI_SUM,
                                                        d_patch_strategy,
                                                        manager);

}

/****************************************************************
* 初始化网格层数据.
*****************************************************************/
void MDLevelIntegrator::initializeLevelData(
    const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int    level_number,
    const double init_data_time,
    const bool   can_be_refined,
    const bool   initial_time,
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
#endif

    NULL_USE(can_be_refined);

    // 初始化数据片.
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
double MDLevelIntegrator::getLevelDt(
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
    const double dt_time,
    const bool initial_time,
    const int  flag_last_dt,
    const double last_dt)
{
    double dt = d_dt_intc->getLevelDt(level,
                                      dt_time,
                                      initial_time,
                                      flag_last_dt,
                                      last_dt,
                                      false);
    return(dt);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int MDLevelIntegrator::advanceLevel(
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
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(predict_dt>=min_dt && predict_dt<=max_dt);
#endif

    NULL_USE(hierarchy_step_number);

    // 显式时间离散格式, 取时间步长等于预测时间步长.
    actual_dt = predict_dt;

    // 调用标准积分函数, 完成1个时间步的积分.
    standardAdvanceLevel(level,
                         current_time,
                         actual_dt,
                         first_step);

    return(1);

}
void MDLevelIntegrator::acceptTimeDependentSolution(
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
    const double new_time,
    const bool deallocate_data)
{
}

/****************************************************************
* 时间步序列中, 积分一个时间步长.
*****************************************************************/
void MDLevelIntegrator::standardAdvanceLevel(
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
    const double current_time,
    const double actual_dt,
    const bool   first_step)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif
    // 规约构件：计算粒子总数目.
    double particle_number = 0;
    d_particle_reduce_intc->reduction( &particle_number, 1, level, current_time, actual_dt, false);
    tbox::plog<< "particle_number="<<particle_number<<endl<<flush;

    // 填充影像区
    d_particle_comm_intc->fillParticles(level, current_time);

    // 更新粒子位置.
    d_particle_comp_intc->computing(level,
                                    current_time,
                                    actual_dt);
    // 迁移粒子.
    d_particle_comm_intc->migrateParticles(level, current_time);
}

/****************************************************************
* 更新网格片中的工作负载.
*****************************************************************/
void MDLevelIntegrator::refreshNonUniformWorkload(
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
    const tbox::Array<double>& measured_load,
    const int    numeric_step,
    const int    workload_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
#endif

    tbox::Pointer< hier::PatchLevel<NDIM> > patch_level = level;
    for ( hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
        tbox::Pointer< hier::Patch<NDIM> > patch = patch_level->getPatch(p());
        d_md_model->computeWorkloadOnPatch(
            *patch, measured_load[0], numeric_step, workload_index);
    }
}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void MDLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nExplicitTimeLevelIntegrator<NDIM>::printClassData..." << endl;
    os << "ExplicitTimeLevelIntegrator<NDIM>: this = "
    << (MDLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;

    os << "d_patch_strategy = "
    << (algs::StandardComponentPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "d_grid_geometry= "
    << (hier::GridGeometry<NDIM>*)d_patch_strategy << endl;

    os << "Integrator Components listed : \n"
    << "  Initialize Integrator Component = "
    << d_init_intc->getName() << " "
    << (algs::IntegratorComponent<NDIM> *)d_init_intc << "\n"
    << "  Step size Integrator Component = "
    << d_dt_intc->getName() << " "
    << (algs::IntegratorComponent<NDIM> *)d_dt_intc << "\n"
    << "  Numerical Integrator Component for particle = "
    << d_particle_comp_intc->getName() << " "
    << (algs::IntegratorComponent<NDIM> *)d_particle_comp_intc << "\n"
    << endl;
}

/****************************************************************
* 从输入文件中读取参数值.
*****************************************************************/
void MDLevelIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{

}

/****************************************************************
* 从重启动数据库中读取参数值.
*****************************************************************/
void MDLevelIntegrator::getFromRestart()
{
}
