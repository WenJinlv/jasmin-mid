//
// 文件名: multiblock_main_Euler_moving_mesh.C
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 170 $
// 修改  : $Date: 2007-06-27 08:44:22 $
// 描述  : 主控程序(多块结构网格拼接, 移动网格方法求解Euler方程组).
//

#include "JASMIN_config.h"

#include "tbox/JASMINManager.h"
#include "tbox/MPI.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/RestartManager.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MemoryUtilities.h"

#include "PatchHierarchy.h"
#include "CartesianCoordinates.h"
#include "MultiblockDeformingGridGeometry.h"
#include "HierarchyTimeIntegrator.h"
#include "JaVisDataWriter.h"
#include "DomainDataWriter.h"

#include "DeformingGridInputUtilities2.h"

#include "MBEulerLevelIntegratorMM.h"
#include "MBEulerMM.h"

using namespace JASMIN;

/***************************************************************************
 * 基于JASMIN框架的多块拼接的结构网格, 采用移动网格r-自适应方法, 解Euler方程组.
 *
 * 主控程序.
 *
 * 该应用实例调用如下的JASMIN框架类:                                        *
 *                                                                          *
 * 核心数据结构类:                                                          *
 *    hier::PatchHierarchy - 以网格层为单位, 管理自适应结构网格.  *
 *                                                                          *
 * 网格几何类:                                                              *
 *     geom::CartesianCoordinates  - 笛卡尔坐标系.                          *
 *     geom::MultiblockDeformingGridGeometry - 移动网格的几何信息.          *
 *     appu::DeformingGridInputUtilities2 - 网格生成辅助工具类
 *                                                                          *
 * 时间积分算法类: algs::HierarchyTimeIntegrator , 需要: *
 *     algs::TimeIntegratorLevelStrategy    - 网格层时间积分算法,需要:      *
 *       algs::IntegratorComponent            - 积分构件.                    *
 *         MBEulerMM                            - 网格片时间积分数值计算子程序,
 **
 *                                             实现各类积分构件.            *
 *
 * 特别说明:
 *    - 该例子采用的移动网格计算方法和Fortran子程序
 *      由北京大学数学学院信息与计算科学系汤华中教授提供.
 *    - 如果对移动网格方法有兴趣, 请参考:
 *         - H.Z.Tang and T.Tang, Adaptive meth methods for one- and
 *           two-dimensional hyperbolic conservation laws, SIAM J.Numer.Anal.,
 *           41(2):487-515, 2003.
 *         - H.Z.Tang, A moving meth method for the Euler flow calculations
 *           using a directional monitor function, COmmunications in
 *           Computational Physics (CiCP), 1(4):656-676, 2006.
 *                                                                          *
 ****************************************************************************
*/

int main(int argc, char* argv[]) {
  /*************************************************************************
   *                                                                       *
   * 一、初始化.                                                           *
   * 包括: JASMIN环境初始化、读入输入参数、                                *
   * 创建应用中需要的JASMIN类对象和自定义类对象、                          *
   * 初始化网格片层次结构和数据片数据等等.                                 *
   *                                                                       *
   *************************************************************************/

  // 1. 初始化MPI和JASMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();

  {
    // 2. 解析命令行参数. 获取计算参数的输入文件名,
    //    重启动目录名及重启动时间步编号.

    string input_filename;
    string restart_read_dirname;
    int restore_num = 0;
    bool is_from_restart = false;
    if ((argc != 2) && (argc != 4)) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                 << "<restart dir> <restore number> [options]\n"
                 << "  options:\n"
                 << "  none at this time" << endl;
      tbox::MPI::abort();
      return (-1);
    } else {
      input_filename = argv[1];
      if (argc == 4) {
        restart_read_dirname = argv[2];
        restore_num = atoi(argv[3]);

        is_from_restart = true;
      }
    }

    // 3. 创建并解析输入文件的计算参数到输入数据库, 称之为根数据库.
    tbox::Pointer<tbox::Database> input_db =
        new tbox::InputDatabase("input_db");
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    // 4. 从根数据库中获得名称为"Main"的子数据库, 并读取控制参数.
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    // 4.1 为日志输出设置参数
    string log_file_name = "euler_mm.log";
    if (main_db->keyExists("log_file_name")) {
      log_file_name = main_db->getString("log_file_name");
    }
    bool log_all_nodes = false;
    if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
    }
    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
    } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
    }
    tbox::plog << "input_filename = " << input_filename << endl;
    tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
    tbox::plog << "restore_num = " << restore_num << endl;

    // 4.2 为可视化输出设置参数
    int viz_dump_interval = 0;
    if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
    }
    string viz_dump_dirname;
    int viz_number_procs_per_file = 1;
    if (viz_dump_interval > 0) {
      if (main_db->keyExists("viz_dump_dirname")) {
        viz_dump_dirname = main_db->getString("viz_dump_dirname");
      }
      if (viz_dump_dirname.empty()) {
        TBOX_ERROR("main(): "
                   << "\nviz_dump_dirname is null ... "
                   << "\nThis must be specified for use with JaVis" << endl);
      }
      if (main_db->keyExists("viz_number_procs_per_file")) {
        viz_number_procs_per_file =
            main_db->getInteger("viz_number_procs_per_file");
      }
    }

    int domain_dump_interval =
        main_db->getIntegerWithDefault("domain_dump_interval", 0);
    string domain_dump_dirname;
    if (domain_dump_interval > 0) {
      if (main_db->keyExists("domain_dump_dirname")) {
        domain_dump_dirname = main_db->getString("domain_dump_dirname");
      }
      if (domain_dump_dirname.empty()) {
        TBOX_ERROR("main(): "
                   << "\ndomain_dump_dirname is null ... "
                   << "\nThis must be specified for use with DomainDataWriter!"
                   << endl);
      }
    }

    // 4.3 为重启动输出设置参数
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

    // 5. 如果是断点续算, 那么打开重启动输入文件.
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    if (is_from_restart) {
      restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                       tbox::MPI::getNodes());
    }

    // 6. 创建计时器.
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    tbox::Pointer<tbox::Timer> t_write_javis_data =
        tbox::TimerManager::getManager()->getTimer(
            "apps::main::write_javis_data");
    tbox::Pointer<tbox::Timer> t_write_restart =
        tbox::TimerManager::getManager()->getTimer(
            "apps::main::write_restart_data");

    // 7. 创建JASMIN框架类的对象, 以及自定义类对象, 并初始化.

    // 7.1 创建笛卡尔坐标系.
    tbox::Pointer<geom::CoordinateSystem<NDIM> > coordinate_system =
        new geom::CartesianCoordinates<NDIM>();

    // 7.2 创建变形结构网格几何.
    tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geometry =
        new geom::MultiblockDeformingGridGeometry<NDIM>(
            "MultiblockGridGeometry", coordinate_system);

    tbox::Pointer<appu::DeformingGridInputUtilities2> grid_tool =
        new appu::DeformingGridInputUtilities2(
            "DeformingGridInputUtilities2",
            input_db->getDatabase("DeformingGridInputUtilities2"),
            grid_geometry);

    // 7.3 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("MultiblockPatchHierarchy",
                                                 grid_geometry);

    // 7.4 创建Euler类对象, 提供单个网格片上求解Euler方程组的数值计算子程序.
    MBEulerMM* euler_model =
        new MBEulerMM("MBEulerMM", input_db->getDatabase("MBEulerMM"),
                      grid_tool, grid_geometry);

    // 7.5 创建网格层时间积分算法.
    MBEulerLevelIntegratorMM* Euler_level_integrator =
        new MBEulerLevelIntegratorMM(
            "MBEulerLevelIntegratorMM",
            input_db->getDatabase("MBEulerLevelIntegratorMM"), euler_model,
            grid_geometry,
            true);  // register to restart manager.

    // 7.6 网格片层次结构显式时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            Euler_level_integrator);

    // 7.7 创建JaVis数据输出器, 并将之注册到Euler类对象.

    tbox::Pointer<tbox::Database> javis_db;
    if (input_db->keyExists("javis")) javis_db = input_db->getDatabase("javis");
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_data_writer =
        new appu::JaVisDataWriter<NDIM>("MBEuler_JaVis_Writer",
                                        viz_dump_dirname,
                                        viz_number_procs_per_file, javis_db);

    tbox::Pointer<appu::DomainDataWriter<NDIM> > domain_writer =
        new appu::DomainDataWriter<NDIM>("MBEuler", domain_dump_dirname);

    // 7.8 注册可视化输出器类对象到用户网格片时间积分算法类.
    euler_model->registerDataWriter(viz_data_writer, domain_writer);

    // 8. *** 初始化网格片层次结构和所有网格片上的数据. ***

    // (1) 创建初始的自适应结构网格, 为数据片调度内存空间,
    //     初始化数据片在初始时刻的值.
    if (!time_integrator->initializeHierarchy()) {
      TBOX_ERROR("\nHierarchy is not successfully initialized. " << endl);
    }

    // (2) 关闭重启动输入文件
    if (is_from_restart) restart_manager->closeRestartFile();

    // (3) 输出初始化信息到日志文件.
    tbox::plog << "\nCheck input data and variables before simulation:" << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

    tbox::plog << "\nCheck Euler data... " << endl;
    euler_model->printClassData(tbox::plog);

    // (4) 输出可视化数据场到可视化输出器.
    t_write_javis_data->start();
    if (viz_dump_interval > 0) {
      viz_data_writer->writePlotData(patch_hierarchy,
                                     time_integrator->getIntegratorStep(),
                                     time_integrator->getIntegratorTime());
    }

    if (domain_dump_interval > 0) {
      domain_writer->writeData(patch_hierarchy,
                               time_integrator->getIntegratorStep(),
                               time_integrator->getIntegratorTime());
    }

    t_write_javis_data->stop();

    /*************************************************************************
     *                                                                       *
     *  二、主循环                                                           *
     *                                                                       *
     *************************************************************************/

    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      int iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::plog << endl
                 << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog << endl
                 << endl;

      // 将数值解推进一个时间步, 返回其时间步长.
      double dt_actual = time_integrator->advanceHierarchy();

      loop_time += dt_actual;

      tbox::plog << endl
                 << endl;
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Dt = " << dt_actual << ", Simulation time is " << loop_time
                 << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog << endl
                 << endl;

      // 在指定时刻, 将数据输出到重启动文件.
      if (write_restart) {
        if ((iteration_num % restart_interval) == 0) {
          t_write_restart->start();
          tbox::RestartManager::getManager()->writeRestartFile(
              restart_write_dirname, iteration_num);
          t_write_restart->stop();
        }
      }

      // 在指定时刻, 将数据输出到绘图文件.
      t_write_javis_data->start();
      if ((viz_dump_interval > 0) && (iteration_num % viz_dump_interval) == 0) {
        viz_data_writer->writePlotData(patch_hierarchy, iteration_num,
                                       loop_time);
      }

      if ((domain_dump_interval > 0) &&
          (iteration_num % domain_dump_interval) == 0) {
        domain_writer->writeData(patch_hierarchy, iteration_num, loop_time);
      }
      t_write_javis_data->stop();
    }

    /*************************************************************************
     *                                                                       *
     *  三、模拟结束. 打印结果，并释放资源.                                  *
     *                                                                       *
     *************************************************************************/

    // 输出计时器统计的时间结果.
    tbox::TimerManager::getManager()->print(tbox::plog);

    tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
    tbox::MemoryUtilities::printMaxMemory(tbox::plog);

    // 释放标准指针.
    if (euler_model) delete euler_model;
    if (Euler_level_integrator) delete Euler_level_integrator;
  }

  // 注销JASMIN系统.
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
