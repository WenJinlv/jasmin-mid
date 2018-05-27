//
// 文件名:    main.C
// 软件包:    JASMIN application
// 版权  :    (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:    $Revision$
// 修改  :    $Date: 2007/03/13 13:10:46 $
// 描述  :    主控程序(求解多群输运方程).
//

#include "JASMIN_config.h"

#include "tbox/JASMINManager.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/RestartManager.h"

#include "CartesianCoordinates.h"
#include "UniRectangularGridGeometry.h"
#include "PatchHierarchy.h"
#include "HierarchyTimeIntegrator.h"
#include "JaVisDataWriter.h"

#include "SnSweepLevelIntegrator.h"
#include "SnSweep.h"

using namespace JASMIN;

/************************************************************************
 *                                                                      *
 * 本示例基于JASMIN框架, 在单层网格上采用Sn方法求解输运方程.            *
 *                                                                      *
 * 特别说明:                                                            *
 *    - 该例子采用的计算方法和Fortran子程序                             *
 *      由北京应用物理与计算数学研究所的阳述林研究员提供.               *
 *                                                                      *
 ***********************************************************************/

/*
 ************************************************************************
 *                                                                      *
 * 运行时, 在命令行中键入输入文件名.                                    *
 *                                                                      *
 *          executable <input file name>                                *
 *                                                                      *
 ************************************************************************
 */

int main(int argc, char* argv[]) {
  // 初始化MPI和JASMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();

  // 解析命令行
  string input_filename;
  string restart_read_dirname;
  int restore_num = 0;
  bool is_from_restart = false;

  if ((argc == 2)) {
    input_filename = argv[1];
  } else if ((argc == 4)) {
    input_filename = argv[1];
    restart_read_dirname = argv[2];
    restore_num = atoi(argv[3]);
    is_from_restart = true;
  } else {
    tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << endl;
    tbox::MPI::abort();
    return (-1);
  }

  if (is_from_restart) {
    tbox::plog << "input_filename = " << input_filename << endl;
    tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
    tbox::plog << "restore_num = " << restore_num << endl;
  }

  // 解析输入文件, 将参数存储到输入数据库中.
  tbox::Pointer<tbox::Database> root_db = new tbox::InputDatabase("root_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, root_db);

  // 从"Main"子数据库中获取日志文件控制参数.
  tbox::Pointer<tbox::Database> main_db = root_db->getDatabase("Main");
  string log_file_name =
      main_db->getStringWithDefault("log_file_name", "snsweep.log");
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
    javis_dump_dirname =
        main_db->getStringWithDefault("javis_dump_dirname", "javis_snsweep");
    javis_number_procs_per_file =
        main_db->getIntegerWithDefault("javis_number_procs_per_file", 1);
  }

  // 设置重启动参数
  int restart_interval = 0;
  if (main_db->keyExists("restart_interval")) {
    restart_interval = main_db->getInteger("restart_interval");
  }

  string restart_write_dirname;
  if (restart_interval > 0) {
    if (main_db->keyExists("restart_write_dirname")) {
      restart_write_dirname = main_db->getString("restart_write_dirname");
    } else {
      TBOX_ERROR("`restart_interval' > 0, but key `restart_write_dirname'"
                 << "not found in input file.");
    }
  }
  const bool write_restart =
      (restart_interval > 0) && !(restart_write_dirname.empty());

  // 如果是断点续算, 那么打开重启动输入文件.
  tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
  if (is_from_restart) {
    restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                     tbox::MPI::getNodes());
  }

  // 从"TimerManager"子数据库获取计时控制参数.
  tbox::TimerManager::createManager(
      root_db->getDatabase("TimerManager"));  // 计时管理器

  {
    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/

    // 笛卡尔坐标系.
    tbox::Pointer<geom::CoordinateSystem<NDIM> > coords =
        new geom::CartesianCoordinates<NDIM>();

    // 均匀矩形网格几何.
    tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geometry =
        new geom::UniRectangularGridGeometry<NDIM>(
            "CartesianGeometry", root_db->getDatabase("CartesianGeometry"),
            coords);

    // 网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

    // 求解输运方程的网格片策略类对象.
    SnSweep* snsweep_model =
        new SnSweep("SnSweep", root_db->getDatabase("SnSweep"), grid_geometry);

    // 求解输运方程的网格层时间积分算法.
    SnSweepLevelIntegrator* snsweep_level_integrator =
        new SnSweepLevelIntegrator(
            "SnSweepLevelIntegrator",
            root_db->getDatabase("SnSweepLevelIntegrator"), snsweep_model,
            grid_geometry);

    // 时间积分算法.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            root_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            snsweep_level_integrator);

    // JaVis数据输出器

    tbox::Pointer<tbox::Database> javis_db;
    if (root_db->keyExists("javis")) javis_db = root_db->getDatabase("javis");
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer;
    if (javis_dump_interval > 0) {
      javis_data_writer = new appu::JaVisDataWriter<NDIM>(
          "JaVis_Writer", javis_dump_dirname, javis_number_procs_per_file,
          javis_db);

      snsweep_model->registerPlotData(javis_data_writer);  //注册待输出的绘图量.
    }

    // 计时器
    tbox::Pointer<tbox::Timer> t_write_javis =
        tbox::TimerManager::getManager()->getTimer("apps::main::write_javis");
    tbox::Pointer<tbox::Timer> t_write_restart =
        tbox::TimerManager::getManager()->getTimer(
            "apps::main::write_restart_data");

    /*******************************************************************************
     *                               初  始  化 *
     *******************************************************************************/

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    //  关闭重启动输入文件
    if (is_from_restart) restart_manager->closeRestartFile();

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

    /************************************************************************************
     *                              时  间  步  进 *
     ************************************************************************************/

    double loop_time = time_integrator->getIntegratorTime();  // 当前时间
    double loop_time_end = time_integrator->getEndTime();     // 终止时间
    int iteration_num = time_integrator->getIntegratorStep();  // 当前时间步数

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      // 将数值解推进一个时间步
      double dt_actual = time_integrator->advanceHierarchy();
      loop_time += dt_actual;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " << iteration_num << endl;
      tbox::pout << "Dt = " << dt_actual << ", Simulation time is " << loop_time
                 << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      iteration_num++;

      // 在指定时刻, 将数据输出到重启动文件.
      if (write_restart) {
        if ((iteration_num % restart_interval) == 0) {
          t_write_restart->start();
          tbox::RestartManager::getManager()->writeRestartFile(
              restart_write_dirname, iteration_num);
          t_write_restart->stop();
        }
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
    if (snsweep_model) delete snsweep_model;
    if (snsweep_level_integrator) delete snsweep_level_integrator;
  }

  // 释放JASMIN和MPI内部资源.
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
