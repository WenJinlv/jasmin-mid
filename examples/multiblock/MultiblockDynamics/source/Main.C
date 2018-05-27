//
// File:       main.C
// 软件包:     JASMIN application
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:    $Revision: 1.13 $
// 修改  :    $Date: 2006/07/19 10:41:31 $
// 描述  :    主控程序(非结构拼接的多块网格, ALE方法求解流体力学方程组).
//

#include "tbox/JASMINManager.h"
#include "tbox/InputManager.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "JaVisDataWriter.h"
#include "CylindricalCoordinates.h"
#include "MultiblockDeformingGridGeometry.h"
#include "PatchHierarchy.h"
#include "HierarchyTimeIntegrator.h"
using namespace JASMIN;

#include <iomanip>
using namespace std;

#include "DynamicsLevelIntegrator.h"
#include "Dynamics.h"

/****************************************************************************
 * 基于JASMIN框架的多块非结构拼接的结构网格,                                *
 * 采用ALE网格r-自适应方法, 求解流体力学方程组.                             *
 *                                                                          *
 * 该应用实例调用如下的JASMIN框架类:                                        *
 *                                                                          *
 * 核心数据结构类:                                                          *
 *    hier::PatchHierarchy - 以网格层为单位, 管理自适应结构网格.  *
 *                                                                          *
 * 网格几何类:                                                              *
 *     geom::CylindricalCoordinates          - 柱坐标系.                    *
 *     geom::MultiblockDeformingGridGeometry - 移动网格的几何信息.          *
 *     appu::DeformingGridInputUtilities2 - 网格生成辅助工具类              *
 *                                                                          *
 * 时间积分算法类:                                                          *
 *     algs::HierarchyTimeIntegrator        - 时间积分算法 *
 *     algs::TimeIntegratorLevelStrategy    - 网格层时间积分算法            *
 *     algs::IntegratorComponent            - 积分构件.                     *
 *     algs::StandardComponentPatchStrategy - 标准构件网格片策略类          *
 *                                                                          *
 ****************************************************************************/
int main(int argc, char* argv[]) {
  // 初始化MPI和JASMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();

  // 解析命令行
  string input_filename;
  string restart_read_dirname;
  int restore_num = 0;
  if ((argc != 2) && (argc != 4)) {
    tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
               << "<restart dir> <restore number> [options]\n"
               << "  options:\n"
               << "  none at this time" << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    // 获取输入文件
    input_filename = argv[1];

    // 参数个数为4时, 获取并打开重启动文件
    if (argc == 4) {
      restart_read_dirname = argv[2];  // 重启动文件的输入路径
      restore_num = atoi(argv[3]);     // 重启动文件的时间序号
      tbox::RestartManager* restart_manager =
          tbox::RestartManager::getManager();  // 重启动管理器
      restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                       tbox::MPI::getNodes());  // 打开文件
    }
  }

  // 解析输入文件, 将参数存储到输入数据库中.
  tbox::Pointer<tbox::Database> root_db = new tbox::InputDatabase("root_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, root_db);

  // 从"Main"子数据库中获取日志文件控制参数.
  tbox::Pointer<tbox::Database> main_db = root_db->getDatabase("Main");
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

  // 从"TimerManager"子数据库获取计时控制参数.
  tbox::TimerManager::createManager(
      root_db->getDatabase("TimerManager"));  // 计时管理器

  {
    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/

    // 柱对称坐标系
    tbox::Pointer<geom::CylindricalCoordinates<NDIM> > cylind_system =
        new geom::CylindricalCoordinates<NDIM>();

    // 多块变形结构网格几何
    tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geometry =
        new geom::MultiblockDeformingGridGeometry<NDIM>(
            "MultiblockDeformingGridGeometry", cylind_system);

    // 变形结构网格几何输入工具
    tbox::Pointer<appu::DeformingGridInputUtilities2> grid_tool =
        new appu::DeformingGridInputUtilities2(
            "DeformingGridInputUtilities2",
            root_db->getDatabase("DeformingGridInputUtilities2"),
            grid_geometry);

    // 网格片层次结构对象.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("DynamicsPatchHierarchy",
                                                 grid_geometry);

    // 求解流体力学方程组的网格片策略类对象.
    Dynamics* dynamics_model = new Dynamics(
        "Dynamics", root_db->getDatabase("Dynamics"), grid_tool, grid_geometry);

    // 求解流体力学方程组的网格层时间积分算法.
    DynamicsLevelIntegrator* dynamics_level_integrator =
        new DynamicsLevelIntegrator(
            "DynamicsLevelIntegrator",
            root_db->getDatabase("DynamicsLevelIntegrator"), dynamics_model,
            grid_geometry);

    // 时间积分算法
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            root_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            dynamics_level_integrator);

    // JaVis数据输出器

    tbox::Pointer<tbox::Database> javis_db;
    if (root_db->keyExists("javis")) javis_db = root_db->getDatabase("javis");
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer;
    if (javis_dump_interval > 0) {
      javis_data_writer = new appu::JaVisDataWriter<NDIM>(
          "JaVis_Writer", javis_dump_dirname, javis_number_procs_per_file,
          javis_db);

      dynamics_model->registerPlotData(
          javis_data_writer);  //注册待输出的绘图量.
    }

    // 计时器
    tbox::Pointer<tbox::Timer> t_write_javis =
        tbox::TimerManager::getManager()->getTimer("apps::main::write_javis");
    tbox::Pointer<tbox::Timer> t_write_restart =
        tbox::TimerManager::getManager()->getTimer("apps::main::write_restart");

    /*******************************************************************************
     *                初 始 化 网 格 片 层 次 结 构 和 物 理 量 *
     *******************************************************************************/

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    // 检查输入数据库参数状态.
    tbox::plog << "\nCheck input data before simulation:" << endl;
    tbox::plog << "Input database..." << endl;
    root_db->printClassData(tbox::plog);

    // 输出初始状态.
    if (javis_dump_interval > 0) {
      t_write_javis->start();
      javis_data_writer->writePlotData(patch_hierarchy,
                                       time_integrator->getIntegratorStep(),
                                       time_integrator->getIntegratorTime());
      t_write_javis->stop();
    }

    // 关闭重启动输入文件.
    tbox::RestartManager::getManager()->closeRestartFile();

    /************************************************************************************
     *                              时  间  步  进 *
     ************************************************************************************/

    double loop_time = time_integrator->getIntegratorTime();  // 当前时间
    double loop_time_end = time_integrator->getEndTime();     // 终止时间
    int iteration_num = time_integrator->getIntegratorStep();  // 当前时间步数

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      tbox::pout << setprecision(12) << endl;
      tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << endl;
      tbox::pout << "Start: Timestep   = " << setw(18) << iteration_num << "\n"
                 << "       Start_Time = " << setw(18) << loop_time << "\n";
      tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << endl;

      // 将数值解推进一个时间步
      double dt_actual = time_integrator->advanceHierarchy();
      loop_time += dt_actual;

      tbox::pout << endl;
      tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << endl;
      tbox::pout << "End: Timestep  = " << setw(18) << iteration_num << "\n"
                 << "     Start_Time= " << setw(18) << loop_time << "\n"
                 << "     Actual_dt = " << setw(18) << dt_actual << "\n"
                 << "     End_Time  = " << setw(18) << loop_time << endl;
      tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << endl;

      iteration_num++;

      // 输出重启动数据.
      if (restart_dump_interval > 0 &&
          (iteration_num % restart_dump_interval) == 0) {
        t_write_restart->start();
        tbox::RestartManager::getManager()->writeRestartFile(
            restart_dump_dirname, iteration_num);
        t_write_restart->stop();
      }

      // 输出可视化数据.
      if ((javis_dump_interval > 0) &&
          (iteration_num % javis_dump_interval) == 0) {
        t_write_javis->start();
        javis_data_writer->writePlotData(patch_hierarchy, iteration_num,
                                         loop_time);
        t_write_javis->stop();
      }
    }

    /************************************************************************************
     *                             结  束  模  拟 *
     ************************************************************************************/

    // 输出计时器统计的时间结果.
    tbox::TimerManager::getManager()->print(tbox::plog);

    // 释放对象.
    if (dynamics_model) delete dynamics_model;
    if (dynamics_level_integrator) delete dynamics_level_integrator;
  }

  // 释放JASMIN和MPI内部资源.
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
