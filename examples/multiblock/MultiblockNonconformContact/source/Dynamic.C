//
// 文件名:  Dynamic.C
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 1.8 $
// 修改  :  $Date: 2007/08/20 01:39:33 $
// 描述  :  用户应用类Dynamic的实现
//

#include "Dynamic.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif

# include <sys/time.h>
# include <unistd.h>

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "GridGeometry.h"
#include "BoxArray.h"
#include "BlockDeformingPatchGeometry.h"
#include "CellData.h"
#include "SideData.h"
#include "NodeData.h"
#include "tbox/IEEE.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "VariableDatabase.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/Tracer.h"
#include "StandardComponentPatchStrategy.h"
#include "tbox/MathUtilities.h"

// 声明Fortran数值子程序的头文件.
#include "DynamicFort.h"

/*
*************************************************************************
* 构造函数, 初始化类Dynamic的数据成员, 创建计算所需变量.                *
* 设置缺省参数后, 如果是从重启动文件恢复模拟, 则该函数调用函数:		*
* getFromRestart()从重启动数据库读取数据.				*
* 最后, 该函数调用函数getFromInput()从输入数据库读取数据.		*
* 对于同一数据, 从输入数据库读取的值将覆盖从重启动数据库读取的值.   	*
*************************************************************************
*/

Dynamic::Dynamic(const string& object_name,
                 tbox::Pointer<tbox::Database> input_db,
                 tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom,
                 tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
                 int module_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!object_name.empty());
  assert(!input_db.isNull());
  assert(!grid_geom.isNull());
  tbox::Tracer tracer("Dynamic:Dynamic");
#endif

  d_object_name = object_name;
  tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

  // 物理模型编号：0: 无效模型.
  //               1: 计算区域包含两个block, 非协调块边界垂直于X轴方向.
  //               2: 计算区域包含两个block, 非协调块边界垂直于Y轴方向.
  //               3: 计算区域包含三个block, 两条非协调块边界都垂直于Y轴方向.
  d_module_id = module_id;

  // 网格生成辅助工具.
  d_grid_tool = grid_tool;

  // 使用指定的网格几何.
  d_grid_geometry = grid_geom;

  d_velocity_in.resizeArray(NDIM);

  // 速度.
  d_velocity   = new pdat::NodeVariable<NDIM, double>("velocity", NDIM);
  // 网格结点坐标.
  d_node_coord = new pdat::NodeVariable<NDIM,double>("node_coord", NDIM);

  /*
   * 从指定的输入数据库/重启动数据库读取数据和控制参数.
   */
  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
  if (is_from_restart)
    {
      getFromRestart();
    }
  getFromInput(input_db, is_from_restart);

  // 向变量数据库注册变量.
  registerModelVariables();
}


/*
*************************************************************************
*                                                                       *
* 类 Dynamic 的析构函数.	 	                                *
*                                                                       *
*************************************************************************
*/

Dynamic::~Dynamic()
{
  tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/*
*************************************************************************
*  									*
* 向变量数据库注册变量, 并获取计算过程中用到的各个数据片的索引号.       *
*									*
*************************************************************************
*/

void Dynamic::registerModelVariables()
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::registerModelVariables");
#endif

  // 获取变量数据库.
  hier::VariableDatabase<NDIM>* variable_db =
    hier::VariableDatabase<NDIM>::getDatabase();

  // 获取当前值上下文, 非协调上下文, 可视化上下文.
  d_current_context    = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
  d_nonconform_context = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NONCONFORM");
  d_plot_context       = d_current_context;

  // 所有变量在CURRENT上下文中的数据片影像区宽度都为0.
  const hier::IntVector<NDIM> zero_ghosts(0);
  // 测试非协调拼接所需数据片的影像区宽度为0.
  const hier::IntVector<NDIM> two_ghosts(0);

  // ***********************************************
  // ************    注册速度   ********************
  // ***********************************************
  // 数据片: (d_velocity, CURRENT), 影像区宽度为zero_ghosts.
  d_velocity_current_id = variable_db->registerVariableAndContext(d_velocity,
                          d_current_context,
                          zero_ghosts);
  // 数据片: (d_velocity, NONCONFORM), 影像区宽度为two_ghosts.
  d_velocity_nonconform_id = variable_db->registerVariableAndContext(d_velocity,
                             d_nonconform_context,
                             two_ghosts);

  // ***************************************************
  // ************    注册网格结点坐标   ****************
  // ***************************************************
  // 数据片: (d_node_coord, CURRENT), 影像区宽度为zero_ghosts.
  d_node_coord_current_id = variable_db->registerVariableAndContext(d_node_coord,
                            d_current_context,
                            zero_ghosts);
  // 数据片: (d_node_coord, NONCONFORM), 影像区宽度为two_ghosts.
  d_node_coord_nonconform_id = variable_db->registerVariableAndContext(d_node_coord,
                               d_nonconform_context,
                               two_ghosts);
}

/*
*************************************************************************
*                                                                       *
* 注册可视化数据片到JaVis输出器, 以采用JaVis对可视化数据文件进行后处理. *
*                                                                       *
*************************************************************************
*/
void Dynamic::registerPlotData(
  tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::registerPlotData");
  assert(!(viz_writer.isNull()));
#endif

  d_javis_writer = viz_writer;

  d_javis_writer->registerNodeCoordinates(d_node_coord_current_id);

#if 0
  hier::VariableDatabase<NDIM>* vardb =
    hier::VariableDatabase<NDIM>::getDatabase();

  d_javis_writer->registerPlotQuantity("Velocity",
                                       "VECTOR",
                                       vardb->mapVariableAndContextToIndex(
                                         d_velocity, d_plot_context));
#endif

}

/********************************************************************
*  初始化积分构件.
*********************************************************************/
void Dynamic::initializeComponent(
  algs::IntegratorComponent<NDIM>* component) const
  {
#ifdef DEBUG_CHECK_ASSERTIONS
    tbox::Tracer tracer("Dynamic::initializeComponent");
    assert(component!=NULL);
#endif

    const string component_name = component->getName();

    if(component_name == "INIT")
      { // 为当前值数据片调度内存.
        component->registerInitPatchData(d_velocity_current_id);
        component->registerInitPatchData(d_node_coord_current_id);
      }
    else if(component_name == "TIME_STEP")
      { // 步长构件, 计算时间步长, 在CURRENT上下文中完成, 不需要通信.
   
      }
    else if(component_name == "UPDATE_NODE_COORDINATE")
      { // 数值构件, 更新网格结点位置坐标, 在CURRENT上下文中完成, 不需要通信.
   
      }
    else if(component_name == "NONCONFORM_CONTACT")
      { // 为非协调拼接构件注册数据片, 这些数据片将从原始网格填充到影像区和非协调影像网格片构成的临时网格层.
        component->registerRefinePatchData(d_node_coord_nonconform_id,
                                           d_node_coord_current_id);
        component->registerRefinePatchData(d_velocity_nonconform_id,
                                           d_velocity_current_id);
        // 注册结点坐标数据片.
        component->registerNodeCoordinatePatchData(d_node_coord_nonconform_id);
        // 为非协调拼接构件注册非协调块边界条件.
        for (int i = 0; i < d_nonconform_boundaries.getSize(); i++)
           component->registerNonconformBoundaryCondition(d_nonconform_boundaries[i]);

        // 注册处理非协调块边界的计算格式的宽度.
        component->registerNonconformOpStencilWidth(d_opstencil_cells);
	
        // 注册警戒网格数目.
        component->registerNonconformGuardWidth(d_guard_cells);
      }
    else if(component_name == "ACCEPT_SOLUTION")
      { // 将数值解从NONCONFORM上下文复制到CURRENT上下文.
        component->registerCopyPatchData(d_velocity_current_id,
                                         d_velocity_nonconform_id);
        component->registerCopyPatchData(d_node_coord_current_id,
                                         d_node_coord_nonconform_id);
      }
    else if(component_name == "NONCONFORM_PATCH_DATA")
      { 
        component->registerPatchData(d_velocity_nonconform_id);
        component->registerPatchData(d_node_coord_nonconform_id);
      }
    else
      {
        TBOX_ERROR(d_object_name << "\n::initializeComponent() : component "
                   << component_name << " is not matched. " << endl);
      }
  }


/*
*************************************************************************
* 当网格片层次结构中增加一个新的网格片时在该网格片内部初始化物理量.	*
* 如果调用该函数的时刻不是初始模拟时刻, 则该函数不做任何操作. 因为在    *
* 其它情况下, 同一网格层进行数据复制就足以初始化一个网格片上的数据.     *
*************************************************************************
*/

void Dynamic::initializePatchData(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const bool initial_time,
                                  const string& component_name)
{
  NULL_USE(time);
#ifdef DEBUG_CHECK_ASSERTIONS

  tbox::Tracer tracer("Dynamic::initializeDataOnPatch");
  assert(component_name == "INIT");
  //tbox::pout << "patch index = " << patch.getPatchNumber() << endl;
#endif

  if (initial_time)
    {
      tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
        patch.getPatchGeometry();

      tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
        patch.getPatchData(d_velocity_current_id);

      tbox::Pointer< pdat::NodeData<NDIM,double> > node_coord =
        patch.getPatchData(d_node_coord_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS

      assert(!velocity.isNull());
      assert(!node_coord.isNull());
#endif

      // 本函数的操作在CURRENT上下文中完成, 所有数据片的影像区宽度均为0.
      const hier::IntVector<NDIM> ghost_cells = velocity->getGhostCellWidth();

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      // 当前多块几何中包含的单块几何数目.
      int sum_blocks = d_grid_geometry->getNumberBlocks();

      // 当前网格块编号.
      int block_id = pgeom->getBlockNumber();

      //tbox::pout << "patch_id = " << patch.getPatchNumber() << endl;

      // 调用Fortran子程序完成初始化操作.
      if (d_module_id == 1) // 第一个物理模型, 非协调块边界垂直于X轴.
        {
          if (sum_blocks != 2)  // 第一个物理模型包含两个块.
            {
              TBOX_ERROR ("module 2 should contain 2 blocks." << endl);
            }
          if (!d_grid_tool.isNull())
            {
              d_grid_tool->generateDeformingMeshForDomain(patch,
                  d_node_coord_current_id);
            }
          else
            {
              TBOX_ERROR ("d_grid_tool is null." << endl);
            }
          initializefield2blockv_(ifirst(0),ilast(0),
                                  ifirst(1),ilast(1),
                                  ghost_cells(0),
                                  ghost_cells(1),
                                  block_id,
                                  d_velocity_in.getPointer(),
                                  velocity->getPointer(),
                                  node_coord->getPointer());
        } 
      else if (d_module_id == 2) // 第二个物理模型, 非协调块边界垂直于Y轴.
        {
          if (sum_blocks != 2)  // 第二个物理模型包含两个块.
            {
              TBOX_ERROR ("module 2 should contain 2 blocks." << endl);
            }
          if (!d_grid_tool.isNull())
            {
              d_grid_tool->generateDeformingMeshForDomain(patch,
                  d_node_coord_current_id);
            }
          else
            {
              TBOX_ERROR ("d_grid_tool is null." << endl);
            }
          initializefield2blockh_(ifirst(0),ilast(0),
                                  ifirst(1),ilast(1),
                                  ghost_cells(0),
                                  ghost_cells(1),
                                  block_id,
                                  d_velocity_in.getPointer(),
                                  velocity->getPointer(),
                                  node_coord->getPointer());
        } 
      else if (d_module_id == 3) // 第三个物理模型, 非协调块边界垂直于Y轴.
        {
          if (sum_blocks != 3)  // 第三个物理模型包含三个块.
            {
              TBOX_ERROR ("module 3 should contain 3 blocks." << endl);
            }
          if (!d_grid_tool.isNull())
            {
              d_grid_tool->generateDeformingMeshForDomain(patch,
                  d_node_coord_current_id);
            }
          else
            {
              TBOX_ERROR ("d_grid_tool is null." << endl);
            }
          initializefield3blockh_(ifirst(0),ilast(0),
                                  ifirst(1),ilast(1),
                                  ghost_cells(0),
                                  ghost_cells(1),
                                  block_id,
                                  d_velocity_in.getPointer(),
                                  velocity->getPointer(),
                                  node_coord->getPointer());
        } 
      else
        {
          TBOX_ERROR ("Wrong d_module_id = " << d_module_id << endl);
        }
    }
}

/*
*************************************************************************
* 在一个网格片上计算并返回时间积分步长.                                 *
*************************************************************************
*/
double Dynamic::getPatchDt(hier::Patch<NDIM>& patch,
                           const double  time,
                           const bool    initial_time,
                           const int     flag_last_dt,
                           const double  last_dt,
                           const string& component_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::getPatchDt");
  assert(component_name == "TIME_STEP");
#endif

  NULL_USE(flag_last_dt);
  NULL_USE(last_dt);

  (void) time;
  (void) initial_time;

  tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
    patch.getPatchGeometry();

  const hier::Index<NDIM> ifirst=patch.getBox().lower();
  const hier::Index<NDIM> ilast =patch.getBox().upper();

  tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
    patch.getPatchData(d_velocity_current_id);

  tbox::Pointer< pdat::NodeData<NDIM,double> > node_coord =
    patch.getPatchData(d_node_coord_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS

  assert(!velocity.isNull());
  assert(!node_coord.isNull());
#endif
  
  return d_dt;
}


/********************************************************************
*  在网格片上更新网格结点位置坐标.
*********************************************************************/
void Dynamic::computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& component_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::computeOnPatch");
#endif

  // 获取当前patch的网格几何.
  tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
    patch.getPatchGeometry();

  // 当前patch所在网格块的序号.
  // int block_id = pgeom->getBlockNumber();
  // int patch_id = patch.getPatchNumber();
  //tbox::pout << "block_id = " << block_id << "  patch_id = " << patch_id << endl;

  const hier::Index<NDIM> ifirst=patch.getBox().lower();
  const hier::Index<NDIM> ilast =patch.getBox().upper();

  tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
    patch.getPatchData(d_velocity_current_id);

  tbox::Pointer< pdat::NodeData<NDIM,double> > node_coord =
    patch.getPatchData(d_node_coord_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS

  assert(!velocity.isNull());
  assert(!node_coord.isNull());
#endif

  const hier::IntVector<NDIM> ghost_cells = velocity->getGhostCellWidth();

  advancecoordinate_(ifirst(0),ilast(0),
                     ifirst(1),ilast(1),
                     ghost_cells(0),
                     ghost_cells(1),
                     dt,
                     velocity->getPointer(),
                     node_coord->getPointer());
}


/*
*************************************************************************
* 找出指定网格片(patch)的所有边界box, 并设置这些边界box中的数据.        *
* 该网格片与物理边界相交的类型(如: 面, 边和结点)由该网格片维护的网格片  *
* 几何类对象获得. 本函数在计算非协调拼接时被调用, 目前主要用于设置      *
* 主非协调拼接外延点信息. 用户                                          *
* 根据当前patch所在的block序号和物理边界类型完成相应操作.               *
*************************************************************************
*/

void Dynamic::setPhysicalBoundaryConditions(
  hier::Patch<NDIM>& patch,
  const double fill_time,
  const hier::IntVector<NDIM>& ghost_width_to_fill,
  const string& component_name)
{
  (void) fill_time;

#ifdef DEBUG_CHECK_ASSERTIONS

  tbox::Tracer tracer("Dynamic::setPhysicalBoundaryConditions");
  assert(component_name == "NONCONFORM_CONTACT");
#endif

  /**
   * 获取当前网格片对应的网格片几何对象.
   */
  const tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
    patch.getPatchGeometry();

  /**
   * 判断当前网格片是否与计算区域边界相接, 如果不相接则返回.
   */
  if (!pgeom->isTouchesPhysicalBoundary())
    return;

  // 获取当前网格片所在的网格块的序号.
  int block_id = pgeom->getBlockNumber();

  /**
   * 获取网格片基本信息.
   */
  const hier::Box<NDIM>& interior = patch.getBox();
  const hier::Index<NDIM>& ifirst = interior.lower();
  const hier::Index<NDIM>& ilast  = interior.upper();

  //tbox::pout << component_name << "   block_id = " << block_id << endl;
  //tbox::pout << component_name << "   ifirst = " << ifirst << endl;
  //tbox::pout << component_name << "   ilast = " << ilast << endl;


  tbox::Pointer< pdat::NodeData<NDIM,double> > node_coordinate =
    patch.getPatchData(d_node_coord_nonconform_id);
  tbox::Pointer< pdat::NodeData<NDIM,double> > node_velocity =
    patch.getPatchData(d_velocity_nonconform_id);

#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!node_coordinate.isNull());
  assert(!node_velocity.isNull());
#endif

  /**
   * 在二维情况下, btype = 1,2, 依次获取"棱"型, "结点"型的boundarybox数组;
   * 在三维情况下, btype = 1,2,3 依次获取"面"型, "棱"型, "结点"型的boundarybox数组.
   */
  const int btype = 1;
  hier::BoxArray<NDIM> fill_boxes;
  tbox::Array<int> location_indexes;
  tbox::Array<int> bdry_conditions;

  hier::IntVector<NDIM> ghost_cells = node_coordinate->getGhostCellWidth();

  pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,
                                      location_indexes,
                                      bdry_conditions,
                                      btype,
                                      interior,
                                      ghost_cells);

  for (int num = 0; num < fill_boxes.size(); num++)
    {
      hier::Box<NDIM> fill_box = fill_boxes(num);
      const hier::Index<NDIM>& ibeg = fill_box.lower();
      const hier::Index<NDIM>& iend = fill_box.upper();
      //tbox::pout << component_name << "   ibeg,iend = " << ibeg << " " << iend << endl;
      //tbox::pout << component_name << "   loc,bdry = " << location_indexes[num]
      //<< " " << bdry_conditions[num] << endl;

      /**
       * 对于"面"类型的物理边界条件:
       * 二维情况下, location_indexes[i]表示fill_boxes位于当前patch的哪个方位, 
       * 具体地: 0: x方向左侧. 
       *         1: x方向右侧. 
       *         2: y方向左侧. 
       *         3: y方向右侧.
       * 三维情况下, 该参数在二维情况基础上, 增加以下值: 
       *         4: z方向左侧.
       *         5: z方向右侧.
       */

    }
}

/********************************************************************
*  在网格片上处理非协调拼接.
*********************************************************************/
void Dynamic::computeNonconformContactOnPatch(hier::Patch<NDIM>& patch,
    const int location,
    hier::Patch<NDIM>& ghost_patch,
    const int ghost_location,
    const int bdry_type,
    const double time,
    const double dt,
    const string& component_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::computeNonconformContactOnPatch");
#endif

  // 获取当前patch的网格几何.
  tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
    patch.getPatchGeometry();

  // 当前patch所在网格块的序号.
  int block_id = pgeom->getBlockNumber();

  // 当前patch的索引范围.
  const hier::Index<NDIM> ifirst=patch.getBox().lower();
  const hier::Index<NDIM> ilast =patch.getBox().upper();

  // 当前patch上的数据片.
  tbox::Pointer< pdat::NodeData<NDIM,double> > node_velocity =
    patch.getPatchData(d_velocity_nonconform_id);

  tbox::Pointer< pdat::NodeData<NDIM,double> > node_coordinate =
    patch.getPatchData(d_node_coord_nonconform_id);

#ifdef DEBUG_CHECK_ASSERTIONS

  assert(!node_velocity.isNull());
  assert(!node_coordinate.isNull());
#endif

  tbox::plog << "block_id = " << block_id << endl;
  tbox::plog << "patch_id = " << patch.getPatchNumber() << endl;
  tbox::plog << "patch box = " << patch.getBox() << endl;
  tbox::plog << "location = " << location << endl;
  tbox::plog << "bdry_type = " << bdry_type << endl;
#if 0
  printPatchData(patch,
                 d_node_coord_nonconform_id,
                 1,  // 数据片类型: 0: 中心量; 1: 结点量.
                 1); // 要打印的变量深度.
#endif
  // 计算新时刻的结点坐标在NONCONFORM上下文中进行.
  const hier::IntVector<NDIM> ghost_cells = node_velocity->getGhostCellWidth();

  // 获取当前patch的网格几何.
  tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > ghost_pgeom =
    ghost_patch.getPatchGeometry();

  // 当前非协调影像patch所在网格块的序号.
  int gblock_id = ghost_pgeom->getBlockNumber();

  // 当前非协调影像patch的索引范围.
  const hier::Index<NDIM> gifirst = ghost_patch.getBox().lower();
  const hier::Index<NDIM> gilast  = ghost_patch.getBox().upper();

  tbox::Pointer< pdat::NodeData<NDIM,double> > ghost_node_velocity =
    ghost_patch.getPatchData(d_velocity_nonconform_id);

  tbox::Pointer< pdat::NodeData<NDIM,double> > ghost_node_coordinate =
    ghost_patch.getPatchData(d_node_coord_nonconform_id);

#ifdef DEBUG_CHECK_ASSERTIONS

  assert(!ghost_node_velocity.isNull());
  assert(!ghost_node_coordinate.isNull());
#endif

  tbox::plog << "gblock_id = " << gblock_id << endl;
  tbox::plog << "gpatch_id = " << ghost_patch.getPatchNumber() << endl;
  tbox::plog << "gpatch box = " << ghost_patch.getBox() << endl;
  tbox::plog << "glocation = " << ghost_location << endl;
#if 0
  printPatchData(ghost_patch,
                 d_node_coord_nonconform_id,
                 1,  // 数据片类型: 0: 中心量; 1: 结点量.
                 0); // 要打印的变量深度.
#endif

}


/*
*************************************************************************
* 将Dynamic类对象的所有数据成员输出到指定的输出流.                        *
*************************************************************************
*/

void Dynamic::printClassData(ostream &os) const
  {

#ifdef DEBUG_CHECK_ASSERTIONS

    tbox::Tracer tracer("Dynamic::printClassData");
#endif

    os << "\nDynamic::printClassData..." << endl;
    os << "Dynamic: this = " << (Dynamic*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_grid_geometry = "
    << (geom::MultiblockDeformingGridGeometry<NDIM>*)d_grid_geometry << endl;

    os << "Parameters for physical problem ..." << endl;

    os << "Problem description and initial data..." << endl;
  }

/*
*************************************************************************
* 从输入数据库读取数据. 需要注意的是从重启动数据库读取的数据值将被从    *
* 输入数据库读取的数据值覆盖.                                           *
*************************************************************************
*/

void Dynamic::getFromInput(
  tbox::Pointer<tbox::Database> db,
  bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::getFromInput");
  assert(!db.isNull());
#endif

  NULL_USE(is_from_restart);

  // 从输入文件读取所有非协调块边界条件序号.
  if (db->keyExists("nonconform_boundaries"))
    {
      d_nonconform_boundaries = db->getIntegerArray("nonconform_boundaries");
    }
  else
    {
      TBOX_ERROR("No key 'nonconform_boundaries' found in input file " << endl);
    }

  if (db->keyExists("velocity"))
    {
      d_velocity_in = db->getDoubleArray("velocity");
    }
  else
    {
      TBOX_ERROR("No key 'velocity' found in input file " << endl);
    }

  if (db->keyExists("time_step"))
    {
      d_dt = db->getDouble("time_step");
    }
  else
    {
      TBOX_ERROR("No key 'time_step' found in input file " << endl);
    }

  if (db->keyExists("opstencil_cells"))
    {
      d_opstencil_cells = db->getDouble("opstencil_cells");
    }
  else
    {
      TBOX_ERROR("No key 'opstencil_cells' found in input file " << endl);
    }

  if (db->keyExists("guard_cells"))
    {
      d_guard_cells = db->getDouble("guard_cells");
    }
  else
    {
      TBOX_ERROR("No key 'guard_cells' found in input file " << endl);
    }
}

/*
*************************************************************************
* 将当前Dynamic类对象的数据成员输出到重启动数据库, 			*
* 或者从重启动数据库读入的子程序.                                       *
*************************************************************************
*/

void Dynamic::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  db->putInteger("Dynamic_version", Dynamic_version);

}

void Dynamic::getFromRestart()
{
  tbox::Pointer<tbox::Database> root_db =
    tbox::RestartManager::getManager()->getRootDatabase();

  tbox::Pointer<tbox::Database> db;
  if ( root_db->isDatabase(d_object_name) )
    {
      db = root_db->getDatabase(d_object_name);
    }
  else
    {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
    }

  int ver = db->getInteger("Dynamic_version");
  if (ver != Dynamic_version)
    {
      TBOX_ERROR(d_object_name << ": "
                 << "Restart file version different than class version." << endl);
    }
}

/*
*************************************************************************
*                                                                       *
* 将参数patch指定网格片上的数据片以参数prec指定的精度输出到参数os指定的 *
* 输出流中. 其中, 实参prec为数据的打印精度(十进制位数), -1表示使用缺省  *
* 缺精度(例如, double类型和complex类型为12位, float类型为6位).          *
*                                                                       *
*************************************************************************
*/

void Dynamic::printPatchData(hier::Patch<NDIM>& patch,
                             const int patchdata_id,
                             const int patchdata_type,
                             const int patchdata_depth,
                             ostream& os,
                             int prec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  tbox::Tracer tracer("Dynamic::printPatchData");
  TBOX_ASSERT ((patchdata_type == 0) || (patchdata_type == 1));
#endif

  tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > pgeom =
    patch.getPatchGeometry();

  int block_id = pgeom->getBlockNumber();

  const hier::Box<NDIM>& interior = patch.getBox();
  const hier::Index<NDIM>& ifirst = interior.lower();
  const hier::Index<NDIM>& ilast  = interior.upper();

  tbox::Pointer< pdat::CellData<NDIM,double> > cell_patch_data;
  tbox::Pointer< pdat::NodeData<NDIM,double> > node_patch_data;

  if (patchdata_type == 0)
    cell_patch_data = patch.getPatchData(patchdata_id);
  else
    node_patch_data = patch.getPatchData(patchdata_id);

#ifdef DEBUG_CHECK_ASSERTIONS
  if (patchdata_type == 0)
    TBOX_ASSERT(cell_patch_data);
  else
    TBOX_ASSERT(node_patch_data);
#endif

  //const int level_id = patch.getPatchLevelNumber();
  //const int patch_id = patch.getPatchNumber();

  //tbox::plog << "level id = " << level_id << endl;
  //tbox::plog << "patch id = " << patch_id << endl;
  //tbox::pout << "block_id = " << block_id << endl;

  //const hier::IntVector<NDIM> ghost_cells = velocity->getGhostCellWidth();

  //  const hier::Box<NDIM>& ghostbox = density->getGhostBox();

  //tbox::plog << "patch range: " << ifirst(0) << " " << ilast(0)
  //           << " " << ifirst(1) << " " << ilast(1) << endl;

  if (patchdata_type == 0)
    cell_patch_data->print(interior,
                           patchdata_depth,
                           os,
                           prec);
  else
    node_patch_data->print(interior,
                           patchdata_depth,
                           os,
                           prec);
}
