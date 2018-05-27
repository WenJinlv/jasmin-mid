//
// 文件名: main_LinAdv.C
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.5 $
// 修改  : $Date: 2007/03/13 13:10:46 $
// 描述  : 主控程序(单层矩形结构网格上, 求解线性对流对流问题).
//

#include "CartesianCoordinates.h"
#include "FederationDeformingGridGeometry.h"
#include "FederationPatchHierarchy.h"
#include "HierarchyTimeIntegrator.h"
#include "JaVisDataWriter.h"
#include "tbox/JASMINManager.h"
#include "tbox/InputManager.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
using namespace JASMIN;

#include "LinAdv.h"
#include "LinAdvLevelIntegrator.h"

#include "LinAdvFederation.h"
#include "LinAdvFederationLevelIntegrator.h"

#include "FederationJaVisDataWriter.h"

/*!
 *************************************************************************
 *
 * @brief 基于JASMIN框架的联邦单层单块均匀矩形结构网格, 求解线性对流问题.
 * 计算区域被分割成左中右三个计算子区域，每个邦员负责计算一个子区域.
 *
 * 该函数分以下几个步骤:
 * -# 预处理: 初始化MPI和JASMIN环境, 解析输入文件, 读取主程序控制参数;
 * -# 创建网格片层次结构和积分算法类对象, 主要包括:
 *    -# 笛卡尔坐标系 geom::CartesianCoordinates<NDIM> ,
 *       联邦均匀矩形网格几何 geom::FederationDeformingGridGeometry<NDIM>,
 *       联邦网格片层次结构(单块) hier::FederationPatchHierarchy<NDIM>
 *    -# 网格片积分算法 LinAdvFederation, LinAdv
 *       - 应用级: 提供求解线性对流方程的数值计算子程序
 *    -# 网格层积分算法 LinAdvFederationLevelIntegrator, LinAdvLevelIntegrator
 *       - 应用级: 提供基于网格层的线性对流方程求解流程
 *    -# 网格片层次结构时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
 * -# 初始化网格片层次结构和物理量数据片;
 * -# 循环: 时间步进;
 * -# 后处理: 释放应用类对象, 释放JASMIN和MPI内部资源.
 *
 ************************************************************************
 */
int main(int argc, char* argv[]) {
  /*******************************************************************************
   *                               预  处  理 *
   *******************************************************************************/
  // 初始化MPI和JASMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();

  // 解析命令行参数:
  bool is_from_restart = false;  // 此次运行是否为重启动.
  string input_filename;         // 输入文件名.
  string restart_read_dirname;   // 重启动文件所在路径
  int restore_num;               // 重启动时间序号
  if ((argc != 2) && (argc != 4)) {
    tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
               << "<restart dir> <restore number> [options]\n"
               << "  options:\n"
               << "  none at this time" << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    input_filename = argv[1];

    if (argc == 4) {  // 重启动情形
      is_from_restart = true;
      restart_read_dirname = argv[2];  // 重启动文件所在路径
      restore_num = atoi(argv[3]);     // 重启动时间序号
    }
  }

  tbox::plog << "input_filename = " << input_filename << endl;
  if (is_from_restart) {
    tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
    tbox::plog << "restore_num = " << restore_num << endl;
  }

  // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
  tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

  // 从根数据库中获得名称为"Main"的子数据库.
  tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

  // 从"Main"子数据库中获取日志文件控制参数.
  string log_file_name = main_db->getString("log_file_name");
  bool log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", false);
  if (log_all_nodes) {
    tbox::PIO::logAllNodes(log_file_name);
  } else {
    tbox::PIO::logOnlyNodeZero(log_file_name);
  }

  // 从"Main"子数据库中获取可视化输出控制参数.
  int javis_dump_interval =
      main_db->getIntegerWithDefault("javis_dump_interval", 0);
  string javis_dump_dirname;
  int javis_number_procs_per_file;
  if (javis_dump_interval > 0) {
    javis_dump_dirname = main_db->getString("javis_dump_dirname");
    javis_number_procs_per_file =
        main_db->getIntegerWithDefault("javis_number_procs_per_file", 1);
  }

  // 从"Main"子数据库中获取重启动输出控制参数.
  int restart_dump_interval =
      main_db->getIntegerWithDefault("restart_dump_interval", 0);
  string restart_dump_dirname;
  if (restart_dump_interval > 0) {
    restart_dump_dirname = main_db->getString("restart_dump_dirname");
  }

  // 打开重启动输入文件
  tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
  if (is_from_restart) {
    restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                     tbox::MPI::getNodes());
  }

  // 创建计时器.
  tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
  tbox::Pointer<tbox::Timer> t_write_javis_data =
      tbox::TimerManager::getManager()->getTimer(
          "apps::main::write_javis_data");
  tbox::Pointer<tbox::Timer> t_write_restart =
      tbox::TimerManager::getManager()->getTimer(
          "apps::main::write_restart_data");

  {
    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/
    // (1) 创建笛卡尔坐标系.
    tbox::Pointer<geom::CoordinateSystem<NDIM> > coords =
        new geom::CartesianCoordinates<NDIM>();

    // (2) 创建联邦均匀矩形网格几何.
    tbox::Pointer<geom::FederationDeformingGridGeometry<NDIM> > grid_geometry =
        new geom::FederationDeformingGridGeometry<NDIM>(
            "FederationCartesianGeometry",
            input_db->getDatabase("FederationCartesianGeometry"), coords);

    // (3) 创建联邦网格片层次结构.
    tbox::Pointer<hier::FederationPatchHierarchy<NDIM> > patch_hierarchy =
        new hier::FederationPatchHierarchy<NDIM>("FederationPatchHierarchy",
                                                 grid_geometry);

    // (4) 创建网格片时间积分算法类.
    const int nfederals = 2;
    LinAdv* linadv_advection_model[nfederals];
    for (int fd = 0; fd < nfederals; fd++) {
      linadv_advection_model[fd] = new LinAdv(
          grid_geometry->getFederalGeometry(fd)->getName() + "_model",
          input_db->getDatabase("LinAdv"));
    }
    LinAdvFederation* linadv_advection_federation_model =
        new LinAdvFederation("LinAdvFederation", grid_geometry);

    // (5) 创建网格层时间积分算法类.
    LinAdvLevelIntegrator* linadv_level_integrator[nfederals];
    for (int fd = 0; fd < nfederals; fd++) {
      linadv_level_integrator[fd] = new LinAdvLevelIntegrator(
          "LinAdvLevel_" + grid_geometry->getFederalGeometry(fd)->getName(),
          linadv_advection_model[fd]);
    }
    LinAdvFederationLevelIntegrator* linadv_level_federation_integrator =
        new LinAdvFederationLevelIntegrator(
            "LinAdvLevelIntegrator", grid_geometry->getNumberFederals(),
            linadv_level_integrator, linadv_advection_federation_model);

    // (6) 创建网格片层次结构时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            linadv_level_federation_integrator);

    // (7) JaVis可视化输出器类.
    tbox::Pointer<tbox::Database> javis_db;
    if (input_db->keyExists("javis")) javis_db = input_db->getDatabase("javis");
    tbox::Pointer<appu::FederationJaVisDataWriter<NDIM> > federal_javis_writer =
        new appu::FederationJaVisDataWriter<NDIM>(
            patch_hierarchy, "federal_javis_datawriter", javis_dump_dirname,
            javis_number_procs_per_file, javis_db);

    if (javis_dump_interval > 0) {
      for (int fd = 0; fd < patch_hierarchy->getNumberFederals(); fd += 1) {
        tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer =
            federal_javis_writer->getDataWriter(fd);
        linadv_advection_model[fd]->registerPlotData(javis_data_writer);
      }
    }

    /*******************************************************************************
     *                初 始 化 网 格 片 层 次 结 构 和 物 理 量 *
     *******************************************************************************/
    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    // 输出初始化信息到日志文件.
    tbox::plog << "\nCheck input data and variables before simulation:" << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

    // 输出可视化数据场到可视化输出器.
    if (javis_dump_interval > 0) {
      t_write_javis_data->start();
      federal_javis_writer->writePlotData(patch_hierarchy,
                                          time_integrator->getIntegratorStep(),
                                          time_integrator->getIntegratorTime());
      t_write_javis_data->stop();
    }

    // 关闭重启动输入文件
    if (is_from_restart) restart_manager->closeRestartFile();

    /************************************************************************************
     *                              时  间  步  进 *
     ************************************************************************************/

    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();

    int iteration_num = time_integrator->getIntegratorStep();

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      // 将数值解推进一个时间步, 返回其时间步长.
      double dt_actual = time_integrator->advanceHierarchy();

      loop_time += dt_actual;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Dt = " << dt_actual << ", Simulation time is " << loop_time
                 << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      // 输出重启动数据.
      if (restart_dump_interval > 0 &&
          (iteration_num % restart_dump_interval) == 0) {
        t_write_restart->start();
        restart_manager->writeRestartFile(restart_dump_dirname, iteration_num);
        t_write_restart->stop();
      }

      // 输出可视化数据.
      if ((javis_dump_interval > 0) &&
          (iteration_num % javis_dump_interval) == 0) {
        t_write_javis_data->start();
        federal_javis_writer->writePlotData(
            patch_hierarchy, time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime());

        t_write_javis_data->stop();
      }
    }

    /************************************************************************************
     *                               模  拟  结  束 *
     ************************************************************************************/
    // 释放对象.
    for (int fd = 0; fd < grid_geometry->getNumberFederals(); fd++) {
      if (linadv_level_integrator[fd]) delete linadv_level_integrator[fd];
      if (linadv_advection_model[fd]) delete linadv_advection_model[fd];
    }

    if (linadv_level_federation_integrator)
      delete linadv_level_federation_integrator;
    if (linadv_advection_federation_model)
      delete linadv_advection_federation_model;
  }

  // 输出计时器统计的时间数据.
  tbox::TimerManager::getManager()->print(tbox::plog);

  // 释放JASMIN和MPI内部资源.
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
