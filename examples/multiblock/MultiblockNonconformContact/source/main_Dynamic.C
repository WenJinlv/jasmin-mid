//
// 文件名: main_Dynamic.C
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.1.1.1 $
// 修改  : $Date: 2007/07/12 13:53:17 $
// 描述  : 演示非协调拼接构件使用方法的主程序.
//

#include "JASMIN_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// JASMIN框架中基本模块的头文件.

#include "tbox/JASMINManager.h"
#include "Box.h"
#include "tbox/Array.h"
#include "BoxList.h"
#include "tbox/Database.h"
#include "Index.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
#include "PatchLevel.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Statistician.h"
#include "VariableDatabase.h"
#include "tbox/MemoryUtilities.h"

// JASMIN中主要算法和数据结构的头文件.

#include "CartesianCoordinates.h"
#include "CylindricalCoordinates.h"
#include "SphericalCoordinates.h"
#include "MultiblockDeformingGridGeometry.h"
#include "HierarchyTimeIntegrator.h"
#include "PatchHierarchy.h"
#include "TimeIntegratorLevelStrategy.h"
#include "JaVisDataWriter.h"

//用户实现的描述应用相关算法和数据结构的头文件.

#include "Dynamic.h"
#include "DynamicLevelIntegrator.h"

using namespace JASMIN;

/*
 * @brief 演示非协调拼接构件使用方法的主程序.
 *
 * 该软件包由以下三部分组成:
 *  - JASMIN框架中的一些类;
 *  - 用户实现的与应用问题相关的类;
 *  - 用户实现的数值子程序.
 *
 * 其中, JASMIN框架提供的类有:
 *
 * hier::PatchHierarchy<NDIM>
 *   - 网格片层次结构, 以网格层为单位, 管理自适应结构网格.
 * geom::CartesianCoordinates<NDIM>
 *   - 提供对笛卡尔坐标系的支持.
 * geom::CylindricalCoordinates<NDIM>
 *   - 提供对柱坐标系的支持.
 * geom::SphericalCoordinates<NDIM>
 *   - 提供对球坐标系的支持.
 * algs::HierarchyTimeIntegrator<NDIM>
 *   - 在网格片层次结构上, 完成时间积分. 并在时间积分过程中.
 *
 * 用户实现的与应用问题相关的类包括以下两个:
 * DynamicLevelIntegrator
 *   - 网格层时间积分算法类.
 * Dynamic
 *   - 网格片时间积分算法类, 在单个网格片上调用Fortran子程序完成计算.
 *     同时创建并管理所有及其数据片.
 *
 * 用户实现的数值子程序见目录: fortran/ .
 *
 ************************************************************************
 */

/*
 **************************************************************************
 * 运行该程序, 必须在命令行提供输入文件名和重启动信息                     *
 *								          *
 * 对于非重启动的情况, 命令如下:				          *
 * 串行执行: 可执行文件 <输入文件名>				          *
 * 并行执行: mpirun -np N 可执行文件 <输入文件名>		          *
 *									  *
 * 对于重启动的情况, 命令如下:					          *
 * 串行执行: 可执行文件 <输入文件名> <重启动文件名> <重启动时间步>        *
 * 并行执行: 								  *
 *    mpirun -np N可执行文件 <输入文件名> <重启动文件名> <重启动时间步>   *
 *								          *
 **************************************************************************
 */

int main(int argc, char* argv[]) {
  /*
   * 初始化MPI和JASMIN环境, 打开日志文件开关并检测命令行提供的参数是否正确.
   */
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();

  string input_filename;
  string restart_read_dirname;
  int restart_num = 0;

  bool is_from_restart = false;

  if ((argc != 2) && (argc != 4)) {
    tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
               << "<restart dir> <restart number> [options]\n"
               << "  options:\n"
               << "  none at this time" << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    input_filename = argv[1];
    if (argc == 4) {
      restart_read_dirname = argv[2];
      restart_num = atoi(argv[3]);

      is_from_restart = true;
    }
  }

  /*
   * 向日志文件中输出输入文件名, 重启动目录以及重启动时间步.
   */
  tbox::plog << "input_filename = " << input_filename << endl;
  tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
  tbox::plog << "restart_num = " << restart_num << endl;

  /*
   * 创建与关键字"input_db"关联的输入数据库,
   * 解析输入数据文件(input_filename)并将其中的数据存储到所创建的数据库中.
   */
  tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

  /*
   * 从数据库(input_db)中获取与关键字"Main"相关联的子数据库,
   * 首先, 从该子数据库读入用于输出日志和可视化数据的控制信息, 如: 日志文件名,
   * 存放可视化文件的目录及文件名等. 其次, 如果在命令行提供了正确的重启动信息,
   * 并且从输入数据库获取的重启动间隔(restart_interval)非零,
   * 则创建一个重启动数据库.
   */

  tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

  string log_file_name = "dynamic.log";
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

  int viz_dump_interval = 0;
  if (main_db->keyExists("viz_dump_interval")) {
    viz_dump_interval = main_db->getInteger("viz_dump_interval");
  }

  tbox::Array<string> viz_writer(1);
  string javis_dump_filename;
  string javis_dump_dirname;

  bool uses_javis = false;
  int javis_number_procs_per_file = 1;

  if (main_db->keyExists("viz_writer")) {
    viz_writer = main_db->getStringArray("viz_writer");
  }

  for (int i = 0; i < viz_writer.getSize(); i++) {
    if (viz_writer[i] == "JaVis") uses_javis = true;
  }

  if (viz_dump_interval > 0) {
    if (main_db->keyExists("viz_dump_filename")) {
      javis_dump_filename = main_db->getString("viz_dump_filename");
      if (uses_javis && javis_dump_filename.empty()) {
        TBOX_ERROR("main(): "
                   << "\njavis_dump_filename is null ... "
                   << "\nThis must be specified for use with JaVis" << endl);
      }
    }
    if (main_db->keyExists("viz_dump_dirname")) {
      javis_dump_dirname = main_db->getString("viz_dump_dirname");
      if (uses_javis && javis_dump_dirname.empty()) {
        TBOX_ERROR("main(): "
                   << "\njavis_dump_dirname is null ... "
                   << "\nThis must be specified for use with JaVis" << endl);
      }
    }
  }
  if (uses_javis && main_db->keyExists("javis_number_procs_per_file")) {
    javis_number_procs_per_file =
        main_db->getInteger("javis_number_procs_per_file");
  }

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

  /*
   * 如果采用多层自适应网格计算, 则从输入文件获取是否采用时间细化模式.
   */
  bool use_refined_timestepping = true;
  if (main_db->keyExists("timestepping")) {
    string timestepping_method = main_db->getString("timestepping");
    if (timestepping_method == "SYNCHRONIZED") {
      use_refined_timestepping = false;
    }
  }

  /*
   * 是否采用隐式格式计算.
   */
  bool use_implicit_discretization = false;
  if (main_db->keyExists("discretization")) {
    string discretization_method = main_db->getString("discretization");
    if (discretization_method == "IMPLICIT") {
      use_implicit_discretization = true;
    }
  }

  const bool write_restart =
      (restart_interval > 0) && !(restart_write_dirname.empty());

  /*
   * 获取重启动管理器和重启动root数据库.
   * 如果是从重启动文件恢复模拟, 则打开重启动文件.
   */

  tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

  if (is_from_restart) {
    restart_manager->openRestartFile(restart_read_dirname, restart_num,
                                     tbox::MPI::getNodes());
  }

  /*
   * 创建计时器管理器, 在程序的运行过程中统计时间性能.
   * 该管理器管理的计时器从输入文件读取. 时间信息将存储在重启动文件中.
   * 如果是从重启动文件恢复模拟, 则计时器将自动被初始化到上一次模拟停止时的状态.
   * 用户也可以调用tbox::TimerManager::resetAllTimers()显式地设置计时器.
   */

  tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

  // 变量module_id表示不同的物理模型.
  int module_id = 0;
  if (main_db->keyExists("module_id")) {
    module_id = main_db->getInteger("module_id");
  } else {
    TBOX_ERROR("No key `module_number' found in input file " << endl);
  }

  tbox::Pointer<geom::CoordinateSystem<NDIM> > coordinate_system;

  // 缺省情况下, 坐标系为笛卡尔坐标系.
  int coordinate_system_tag = 0;
  if (main_db->keyExists("coordinate_system_tag")) {
    coordinate_system_tag = main_db->getInteger("coordinate_system_tag");
    if (coordinate_system_tag == 0) {
      // 创建笛卡尔坐标系.
      coordinate_system = new geom::CartesianCoordinates<NDIM>;
    } else if (coordinate_system_tag == 1) {
      // 创建柱坐标系.
      coordinate_system = new geom::CylindricalCoordinates<NDIM>;
    } else if (coordinate_system_tag == 2) {
      // 创建球坐标系.
      coordinate_system = new geom::SphericalCoordinates<NDIM>;
    } else {
      TBOX_ERROR("Wrong coordinate_system_tag: " << coordinate_system_tag
                                                 << endl);
    }
  }

  /**************************************************************
   *   下面创建与应用问题相关的主要算法和数据对象,              *
   *   并从输入文件或者重启动文件或者二者的结合初始化每个对象.  *
   **************************************************************/

  // 多块变形网格几何.
  tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geometry =
      NULL;

  // 变形结构网格几何输入工具.
  tbox::Pointer<appu::DeformingGridInputUtilities2> grid_tool = NULL;

  // 创建变形网格几何.
  grid_geometry = new geom::MultiblockDeformingGridGeometry<NDIM>(
      "MultiblockDeformingGridGeometry", coordinate_system);

  // 创建变形结构网格几何输入工具.
  grid_tool = new appu::DeformingGridInputUtilities2(
      "DeformingGridInputUtilities2",
      input_db->getDatabase("DeformingGridInputUtilities2"), grid_geometry);

  // 创建网格片层次结构.
  tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
      new hier::PatchHierarchy<NDIM>("MultiblockPatchHierarchy",
                                               grid_geometry);

  /**
   * 创建网格片时间积分算法类对象, 该对象在单个网格片上完成计算,
   * 如: 初始化网格片数据, 计算时间步长等. 另外, 创建并管理变量及其数据片.
   */
  Dynamic* dynamic_model =
      new Dynamic("Dynamic", input_db->getDatabase("Dynamic"), grid_geometry,
                  grid_tool, module_id);

  /**
   * 网格层时间积分算法类.
   */
  DynamicLevelIntegrator* dynamic_level_integrator = new DynamicLevelIntegrator(
      "DynamicLevelIntegrator", input_db->getDatabase("DynamicLevelIntegrator"),
      dynamic_model, grid_geometry, patch_hierarchy,
      true,  // register to restart manager.
      use_refined_timestepping, use_implicit_discretization);

  // 网格片层次结构时间积分算法类.
  tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
      new algs::HierarchyTimeIntegrator<NDIM>(
          "HierarchyTimeIntegrator",
          input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
          dynamic_level_integrator);

  // 创建JaVis可视化输出器类对象.

  tbox::Pointer<tbox::Database> javis_db;
  if (input_db->keyExists("javis")) javis_db = input_db->getDatabase("javis");
  tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer =
      new appu::JaVisDataWriter<NDIM>("Dynamic_JaVis_Writer",
                                      javis_dump_dirname,
                                      javis_number_procs_per_file, javis_db);

  // 注册可视化输出器类对象到用户网格片时间积分算法类.
  dynamic_model->registerPlotData(javis_data_writer);

  /******************************************************
   *  创建初始网格, 为数据片调度内存空间, 初始化数据片. *
   ******************************************************/
  if (!time_integrator->initializeHierarchy()) {
    TBOX_ERROR("\nHierarchy is not successfully initialized. " << endl);
  }

  // 打印网格片层次结构信息.
  // patch_hierarchy->recursivePrint(tbox::pout,"",1);

  // 关闭重启动数据库.
  tbox::RestartManager::getManager()->closeRestartFile();

  /*
   * 创建计时器用以测量I/O性能.
   */

  tbox::Pointer<tbox::Timer> t_write_viz =
      tbox::TimerManager::getManager()->getTimer("apps::main::write_viz");
  tbox::Pointer<tbox::Timer> t_write_restart =
      tbox::TimerManager::getManager()->getTimer("apps::main::write_restart");

  // 输出可视化数据场到可视化输出器.
  t_write_viz->start();
  if (viz_dump_interval > 0) {
    javis_data_writer->writePlotData(patch_hierarchy,
                                     time_integrator->getIntegratorStep(),
                                     time_integrator->getIntegratorTime());
  }
  t_write_viz->stop();

  /*
   * 获取当前积分步的时间和积分终止时间.
   */
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

#if 1
    // 每隔一定的时间步间隔, 打印网格片层次结构信息及时间统计结果.
    if ((iteration_num % 10) == 0) {
      // patch_hierarchy->recursivePrint(tbox::pout,"",1);
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::TimerManager::getManager()->print(tbox::plog);
    }
#endif

    /*
     * 在指定的时刻, 将数据输出到重启动文件.
     */
    if (write_restart) {
      if ((iteration_num % restart_interval) == 0) {
        t_write_restart->start();
        tbox::RestartManager::getManager()->writeRestartFile(
            restart_write_dirname, iteration_num);
        t_write_restart->stop();
      }
    }

    /*
     * 在指定的时刻, 将数据输出到.
     */
    t_write_viz->start();
    if ((viz_dump_interval > 0) && (iteration_num % viz_dump_interval) == 0) {
      if (uses_javis) {
        javis_data_writer->writePlotData(patch_hierarchy, iteration_num,
                                         loop_time);
      }
    }
    t_write_viz->stop();
  }

  /*
   * 输出计时器统计出的各种时间开销.
   */
  tbox::TimerManager::getManager()->print(tbox::plog);
#if 0
  /*
   * 输出内存使用情况.
   */
  tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
  tbox::MemoryUtilities::printMaxMemory(tbox::plog);
#endif

  /*
   * 模拟结束后, 释放所有对象.
   */
  patch_hierarchy.setNull();
  grid_geometry.setNull();
  time_integrator.setNull();

  javis_data_writer.setNull();

  if (dynamic_model) {
    delete dynamic_model;
  }

  if (dynamic_level_integrator) {
    delete dynamic_level_integrator;
  }

  /*
   * 结束MPI和JASMIN环境.
   */
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
