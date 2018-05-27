//
// 文件  :	Dynamics.C
// 软件包:	JASMIN application
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 1.29 $
// 修改  :	$Date: 2006/08/08 01:05:45 $
// 描述  :	求解流体力学方程的用户网格片子程序的具体实现.
//

#include "Dynamics.h"
#include "DynamicsFort.h"

#include "BlockDeformingPatchGeometry.h"
#include "CellData.h"
#include "SideData.h"
#include "NodeData.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "PatchNodeDataOpsReal.h"
#include "tbox/IEEE.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <assert.h>

// 表示边界条件的常量
#include "ConstantDefines.h"

/*
*****************************************************************
* 构造函数.
*
* 该函数按如下步骤执行:  
* (1) 从重启动数据库和输入数据库中读取参数, 初始化私有数据成员.
* (2) 注册变量.
*****************************************************************
*/

Dynamics::Dynamics(const string& object_name,
             tbox::Pointer<tbox::Database> input_db,
             tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
             tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
   assert(!grid_tool.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometry = grid_geom;
   d_grid_tool     = grid_tool;
  
   /*
    * 从指定的输入数据库/重启动数据库初始化数据成员.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   registerModelVariables();
}


/*
*************************************************************************
*                                                                       *
* 类 Dynamics 的空析构函数.	 	                                *
*                                                                       *
*************************************************************************
*/

Dynamics::~Dynamics() 
{
}


/*
*************************************************************************
*  									*
* 注册变量.
*
*************************************************************************
*/

void Dynamics::registerModelVariables()
{
   // 为影像区的宽度赋值.
   hier::IntVector<NDIM> zeroghosts(0);
   hier::IntVector<NDIM> oneghosts(1);

   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   /*
    * 创建流体力学的基本物理量.
    */
   d_coordinates = new pdat::NodeVariable<NDIM,double>("coordinates",NDIM);
   d_velocity    = new pdat::NodeVariable<NDIM,double>("velocity"  ,NDIM);
   d_density     = new pdat::CellVariable<NDIM,double>("density" ,1);
   d_energy      = new pdat::CellVariable<NDIM,double>("energy",1);
   d_gamma       = new pdat::CellVariable<NDIM,double>("gamma",1);

   /*
    * 创建流体力学的导出物理量.
    */
   d_speedup     = new pdat::NodeVariable<NDIM,double>("speedup",NDIM);
   d_pressure    = new pdat::CellVariable<NDIM,double>("pressure",1);
   d_mass        = new pdat::CellVariable<NDIM,double>("mass",1);
   d_viscosity   = new pdat::CellVariable<NDIM,double>("viscosity",1);
   d_paq         = new pdat::CellVariable<NDIM,double>("paq",1);

   /*
    * 创建辅助物理量.
    */
   d_side_pressure  = new pdat::SideVariable<NDIM,double>("side_pressure",1);


   // 当前值上下文, 新值上下文, 演算上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");

   // 存储当前时刻的变量值: (coordnates,current), 影像区宽度为0.
   d_coordinates_current_id =
         variable_db->registerVariableAndContext(d_coordinates,
                                                 d_current,
                                                 zeroghosts);

   // 存储新时刻的变量值: (coordinates,new), 影像区宽度为0.
   d_coordinates_new_id =
         variable_db->registerVariableAndContext(d_coordinates,
                                                 d_new,
                                                 zeroghosts);

   // 存储变量的演算值: (coordinatess,scratch), 影像区宽度等于格式宽度.
   d_coordinates_scratch_id =
         variable_db->registerVariableAndContext(d_coordinates,
                                                 d_scratch,
                                                 oneghosts);

   // 存储当前时刻的变量值: (density,current), 影像区宽度为0.
   d_density_current_id =
         variable_db->registerVariableAndContext(d_density,
                                                 d_current,
                                                 zeroghosts);


   // 存储新时刻的变量值: (density,new), 影像区宽度为0.
   d_density_new_id =
         variable_db->registerVariableAndContext(d_density,
                                                 d_new,
                                                 zeroghosts);

   // 存储变量的演算值: (densitys,scratch), 影像区宽度等于格式宽度.
   d_density_scratch_id =
         variable_db->registerVariableAndContext(d_density,
                                                 d_scratch,
                                                 oneghosts);

   // 存储当前时刻的变量值: (velocity,current), 影像区宽度为0.
   d_velocity_current_id =
         variable_db->registerVariableAndContext(d_velocity,
                                                 d_current,
                                                 zeroghosts);

   // 存储新时刻的变量值: (velocity,new), 影像区宽度为0.
   d_velocity_new_id =
         variable_db->registerVariableAndContext(d_velocity,
                                                 d_new,
                                                 zeroghosts);

   // 存储变量的演算值: (velocity,scratch), 影像区宽度等于格式宽度.
   d_velocity_scratch_id =
         variable_db->registerVariableAndContext(d_velocity,
                                                 d_scratch,
                                                 oneghosts);

   // 存储当前时刻的变量值: (energy,current), 影像区宽度为0.
   d_energy_current_id =
         variable_db->registerVariableAndContext(d_energy,
                                                 d_current,
                                                 zeroghosts);


   // 存储新时刻的变量值: (energy,new), 影像区宽度为0.
   d_energy_new_id =
         variable_db->registerVariableAndContext(d_energy,
                                                 d_new,
                                                 zeroghosts);

   d_gamma_current_id =
         variable_db->registerVariableAndContext(d_gamma,
                                                 d_current,
                                                 zeroghosts);

   // 存储变量的演算值: (pressure,scratch), 影像区宽度等于格式宽度.
   d_pressure_scratch_id =
         variable_db->registerVariableAndContext(d_pressure,
                                                 d_scratch,
                                                 oneghosts);
   
   // 存储变量的演算值: (side_pressure,scratch), 影像区宽度等于格式宽度.
   d_side_pressure_scratch_id =
         variable_db->registerVariableAndContext(d_side_pressure,
                                                 d_scratch,
                                                 oneghosts);
   
    // 存储当前时刻的变量值: (viscosity,current), 影像区宽度为1.
   d_viscosity_scratch_id =
         variable_db->registerVariableAndContext(d_viscosity,
                                                 d_current,
                                                 oneghosts);

    // 存储当前时刻的变量值: (paq,scratch), 影像区宽度为1.
   d_paq_scratch_id =
         variable_db->registerVariableAndContext(d_paq,
                                                 d_scratch,
                                                 oneghosts);

     // 存储当前时刻的变量值: (mass,scratch), 影像区宽度为1.
   d_mass_scratch_id =
         variable_db->registerVariableAndContext(d_mass,
					   d_scratch,
					   oneghosts);

     // 存储当前时刻的变量值: (speedup,current), 影像区宽度为0.
   d_speedup_scratch_id =
         variable_db->registerVariableAndContext(d_speedup,
					   d_scratch,
					   oneghosts);

 
}

/*
*************************************************************************
*
* 初始化积分构件.
*
*************************************************************************
*/
void Dynamics::initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件 : 为网格坐标和数值解赋初值.
        intc->registerInitPatchData(d_coordinates_current_id);
        intc->registerInitPatchData(d_velocity_current_id);
        intc->registerInitPatchData(d_density_current_id);
        intc->registerInitPatchData(d_energy_current_id);
        intc->registerInitPatchData(d_gamma_current_id);
   }else if(intc_name=="TIME_STEP_SIZE") { //步长构件: 求时间步长.

   }else if(intc_name=="PRE") { //数值构件: 前处理.

   }else if(intc_name=="SIDE_PRESSURE") { //数值构件: 边压力.
        intc->registerRefinePatchData(d_coordinates_scratch_id,
                                      d_coordinates_current_id);
        intc->registerRefinePatchData(d_density_scratch_id,
                                      d_density_current_id);
        intc->registerRefinePatchData(d_velocity_scratch_id,
                                      d_velocity_current_id);
        intc->registerRefinePatchData(d_mass_scratch_id,
                                      d_mass_scratch_id);
        intc->registerRefinePatchData(d_paq_scratch_id,
                                      d_paq_scratch_id);
       
   }else if(intc_name=="VELOCITY") {//数值构件：计算 
        intc->registerRefinePatchData(d_side_pressure_scratch_id,
                                      d_side_pressure_scratch_id);

   }else if(intc_name=="ENERGY") {//数值构件：计算 

   }else if(intc_name=="ROUNDOFF") {//舍入误差构件
        intc->registerPatchData(d_velocity_new_id);
  
   }else if(intc_name=="ACCEPT_SOLUTION") { //复制构件: 接收数值解.
        intc->registerCopyPatchData(d_coordinates_current_id,
                                    d_coordinates_new_id);
        intc->registerCopyPatchData(d_density_current_id,
                                    d_density_new_id);
        intc->registerCopyPatchData(d_velocity_current_id,
                                    d_velocity_new_id);
        intc->registerCopyPatchData(d_energy_current_id,
                                    d_energy_new_id);

   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 内存构件: 为新值数据片调度内存空间.
        intc->registerPatchData(d_coordinates_new_id);
        intc->registerPatchData(d_density_new_id);
        intc->registerPatchData(d_velocity_new_id);
        intc->registerPatchData(d_energy_new_id);
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 内存构件: 为演算数据片调度内存空间.
        intc->registerPatchData(d_coordinates_scratch_id);
        intc->registerPatchData(d_density_scratch_id);
        intc->registerPatchData(d_velocity_scratch_id);
        intc->registerPatchData(d_pressure_scratch_id);
        intc->registerPatchData(d_speedup_scratch_id);
        intc->registerPatchData(d_paq_scratch_id);
        intc->registerPatchData(d_viscosity_scratch_id);
        intc->registerPatchData(d_mass_scratch_id);
        intc->registerPatchData(d_side_pressure_scratch_id);
   }else{
        TBOX_ERROR("\n::initializeComponent() : component "
                   << intc_name <<" is not matched. "<<endl);
   }

}
 

/*
*************************************************************************
* 在模拟的初始时刻,根据属性区的性质,初始化各个网格片内部的网格坐标,速度,
* 密度和温度.
*
* 该函数仅在初始时刻被调用.在其余情形,可复制旧网格层的变量值完成初始化.
*************************************************************************
*/

void Dynamics::initializePatchData(
             hier::Patch<NDIM>& patch,
             const double data_time,
             const bool initial_time,
             const string& intc_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
TBOX_ASSERT(intc_name=="INIT");
#endif

   (void) data_time;

   if (initial_time) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates_current =
                               patch.getPatchData(d_coordinates_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > density_current =
                               patch.getPatchData(d_density_current_id);
      tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_current =
                               patch.getPatchData(d_velocity_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > energy_current =
                               patch.getPatchData(d_energy_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > gamma_current =
                               patch.getPatchData(d_gamma_current_id);
         
      // 生成初始网格.
      if(!d_grid_tool.isNull()) {
         d_grid_tool->generateDeformingMeshForDomain(patch,
                                                     d_coordinates_current_id);
      }

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      // 初始状态.
      density_current->fill(d_initial_density);
      energy_current ->fill(d_initial_energy);
      gamma_current  ->fill(d_initial_gamma);
          
      for (int d=0;d<NDIM;d++) {
          velocity_current->fill(d_initial_velocity[d],d);
      }  

      const hier::Index<NDIM> center(0);
      tbox::Pointer< hier::BlockPatchGeometry<NDIM> > bgeom =
                                          patch.getPatchGeometry();

      const hier::IntVector<NDIM> ghost_cells(1);
        
      //For Sedov Polar mesh problem

      if( d_model == "SEDOV" ){
         if ( bgeom->getBlockNumber()==0 && ifirst[1] == 0 ) { 
             hier::Index<NDIM> centerlast =ilast;
	     centerlast[1] = 0;
             const hier::Box<NDIM> centerbox(ifirst,centerlast);
             energy_current->fill(d_initial_center_energy, centerbox);
         }
      }else{
         TBOX_ERROR("wrong model name!");
      }
   } 
}


/*
*************************************************************************
* 计算时间步长.
*************************************************************************
*/
double Dynamics::getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name)
{
   const tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates = 
      patch.getPatchData(d_coordinates_current_id);
   const tbox::Pointer< pdat::NodeData<NDIM,double> > velocity = 
      patch.getPatchData(d_velocity_current_id);

   const tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density_current_id);
   const tbox::Pointer< pdat::CellData<NDIM,double> > energy   = 
      patch.getPatchData(d_energy_current_id);

   const hier::IntVector<NDIM> ghost_cells  =  density->getGhostCellWidth();
   const hier::IntVector<NDIM> ifirst_cells = (density->getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (density->getBox()).upper();

   double stabdt = 0.;
   stable_dt_(
             ifirst_cells(0),ilast_cells(0),
             ifirst_cells(1),ilast_cells(1),
             ghost_cells(0),ghost_cells(1),
             coordinates->getPointer(),
             velocity->getPointer(),
             density->getPointer(),
             energy->getPointer(),
             stabdt);

   return(stabdt);
}

/*
*************************************************************************
* 为数值构件完成网格片数值计算.
*************************************************************************
*/
void Dynamics::computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name)
{
   tbox::Pointer< hier::BlockPatchGeometry<NDIM> > patch_geometry =
                                          patch.getPatchGeometry();
   if(intc_name=="PRE") {
      preprocessStateOnPatch(patch,time,dt);
   }else if(intc_name=="SIDE_PRESSURE") {
      computeSidePressureOnPatch(patch,time,dt);
   }else if(intc_name=="VELOCITY") {
      computeSpeedupOnPatch(patch);
      ANCOnPatch(patch);
      computeVelocityOnPatch(patch,time,dt);
   }else if(intc_name=="ENERGY") {
      computeGridDensityEnergyOnPatch(patch,time,dt);
   }else{
        TBOX_ERROR("\n::initializeComponent() : component "
                   << intc_name <<" is not matched. "<<endl);
   }
}

/*
*************************************************************************
* 预处理流体的基本量,包括压力和黏性等.
*************************************************************************
*/
void Dynamics::preprocessStateOnPatch(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const double dt)
{
   (void) time;
   (void) dt;

   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates =
      patch.getPatchData(d_coordinates_current_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity_current_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_scratch =
      patch.getPatchData(d_velocity_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > viscosity =
      patch.getPatchData(d_viscosity_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > density =
      patch.getPatchData(d_density_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > mass =
      patch.getPatchData(d_mass_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > energy =
      patch.getPatchData(d_energy_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > gamma =
      patch.getPatchData(d_gamma_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
      patch.getPatchData(d_pressure_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > paq =
      patch.getPatchData(d_paq_scratch_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(pressure->getGhostCellWidth()==ghost_cells);
   assert(paq->getGhostCellWidth()==ghost_cells);
   assert(viscosity->getGhostCellWidth()==ghost_cells);
#endif 

   
   preprocess_dynamics_state_( 
         ifirst_cells(0),ilast_cells(0),
         ifirst_cells(1),ilast_cells(1),
         ghost_cells(0),ghost_cells(1),
         coordinates -> getPointer(),    // 结点坐标
         density     -> getPointer(),    // 网格密度     
         mass        -> getPointer(),    // 网格质量  
         velocity    -> getPointer(),    // 结点速度 
         energy      -> getPointer(),    // 能量 
         viscosity   -> getPointer(),    // 人为粘力
         paq         -> getPointer(),    // 压力+人为粘力
         pressure    -> getPointer(),    // 压力
         gamma       -> getPointer(),
	 d_a0,
	 d_b0         ); 
}

/*
*************************************************************************
* 计算边压力.
*************************************************************************
*/

void Dynamics::computeSidePressureOnPatch(hier::Patch<NDIM>& patch,
                                 const double time, 
                                 const double dt)
{

   tbox::Pointer< pdat::CellData<NDIM,double> > mass  = 
      patch.getPatchData(d_mass_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > paq = 
      patch.getPatchData(d_paq_scratch_id);

   tbox::Pointer< pdat::SideData<NDIM,double> > side_pressure = 
      patch.getPatchData(d_side_pressure_scratch_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!mass.isNull());
   assert(paq->getGhostCellWidth()==ghost_cells);
   assert(mass->getGhostCellWidth()==ghost_cells);
#endif

   side_pressure->fillAll(0.0);
   // 调用Fortran子程序,计算边压力.
   compute_side_pressure_(
             ifirst_cells(0),ilast_cells(0),
             ifirst_cells(1),ilast_cells(1),
             ghost_cells(0),ghost_cells(1),
             mass   ->getPointer(),        // 质量
             paq    ->getPointer(),        // 压力
             side_pressure->getPointer(0),  // 边压力
             side_pressure->getPointer(1)   // 边压力
             );
}


/*
*************************************************************************
* 计算网格片内部结点的加速度.
*************************************************************************
*/

void Dynamics::computeSpeedupOnPatch(hier::Patch<NDIM>& patch)
{

    tbox::Pointer< pdat::NodeData<NDIM,double> > speedup
                         = patch.getPatchData(d_speedup_scratch_id);
    speedup->fillAll(0.0);

    math::PatchNodeDataOpsReal<NDIM,double> ops;
    const hier::Box<NDIM> &patch_box = patch.getBox();
    hier::IntVector<NDIM> oneghost(1);
    hier::Box<NDIM> ghost_box = hier::Box<NDIM>::grow(patch_box,oneghost);

    const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom =
	                patch.getPatchGeometry();
    if ( !pgeom->isTouchesBlockSingularBoundary() ) { // 不接触奇异点

	computeSpeedupOnBox(patch,ghost_box);
        double alpha = 0.25;
        ops.scale(speedup, alpha, speedup, patch_box );

    } else {
        // 获取奇异点信息
        hier::BoxArray<NDIM> fill_boxes;
        tbox::Array<int> location_indexes;
        tbox::Array<tbox::Array<tbox::Pointer<hier::Patch<NDIM > > > > singu_patches;
        tbox::Array< tbox::Array< int > > blocks;
	int codim = NDIM;
        const hier::IntVector<NDIM> gcw(1);
	pgeom->getBlockInternalSingularBoundaryFillBoxes (
                                           fill_boxes, 
		                           location_indexes, 
					   singu_patches, 
					   blocks, 
					   codim, 
					   patch_box, 
					   gcw);

        // 处理非奇异点
        hier::BoxList<NDIM> box_list(ghost_box);
        box_list.removeIntersections(fill_boxes);
        hier::BoxArray<NDIM> boxes(box_list);
        for ( int i=0; i<boxes.size(); i++ ) {
            computeSpeedupOnBox(patch,boxes(i));
        }

        double alpha = 0.25;
        ops.scale(speedup, alpha, speedup, patch_box );

        // 处理奇异点
        for ( int i=0; i<fill_boxes.size(); i++ ) {	
           speedup->fillAll(0.0, fill_boxes(i));
	   for(int n=0; n<singu_patches[i].size(); n++) {
	      tbox::Pointer<hier::Patch<NDIM> > singu_patch = singu_patches[i][n];
	      tbox::Pointer< pdat::NodeData<NDIM,double> > speedup_sing 
                         = singu_patch->getPatchData(d_speedup_scratch_id);
	      speedup_sing->fillAll(0.0);
	      computeSpeedupOnBox(*singu_patch,
		                  singu_patch->getBox());
              ops.add(speedup,speedup,speedup_sing,fill_boxes(i));
	   }

           double alpha = 1.0 / singu_patches[i].size();
           ops.scale(speedup, alpha, speedup, fill_boxes(i) );
        }
    }

}


void Dynamics::computeSpeedupOnBox(hier::Patch<NDIM>& patch,
                                 const hier::Box<NDIM>& box)
{

   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates =
      patch.getPatchData(d_coordinates_scratch_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > speedup =
      patch.getPatchData(d_speedup_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > paq = 
      patch.getPatchData(d_paq_scratch_id);
   tbox::Pointer< pdat::SideData<NDIM,double> > side_pressure = 
      patch.getPatchData(d_side_pressure_scratch_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coordinates.isNull());
   assert(!speedup.isNull());
   assert(!density.isNull());
   assert(paq->getGhostCellWidth()==ghost_cells);
#endif


   compute_speedup_iga_(
             ifirst_cells(0),ilast_cells(0),
             ifirst_cells(1),ilast_cells(1),
             ghost_cells(0),ghost_cells(1),
	     box.lower(0), box.upper(0),
	     box.lower(1), box.upper(1),
             coordinates->getPointer(),    // 结点坐标
             density->getPointer(),        // 密度
             paq    ->getPointer(),       // 压力
             side_pressure    ->getPointer(0),       // 压力
             side_pressure    ->getPointer(1),       // 压力
             speedup->getPointer());        // 结点加速度



}



/*
*************************************************************************
* ANC方法修正速度.   				                        *
*************************************************************************
*/

void Dynamics::ANCOnPatch(hier::Patch<NDIM>& patch )
{

    if(d_anc < 1.0e-6) return;

    tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_scratch
       = patch.getPatchData(d_velocity_scratch_id);
    tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_current
       = patch.getPatchData(d_velocity_current_id);

    const hier::Box<NDIM> &patch_box = patch.getBox();
    tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_tmp 
	= new pdat::NodeData<NDIM,double>(patch_box,NDIM,1); //保存修正量的临时数据片
    velocity_tmp->fillAll(0.0);

    hier::IntVector<NDIM> oneghost(1);
    hier::Box<NDIM> ghost_box = hier::Box<NDIM>::grow(patch_box,oneghost);
    math::PatchNodeDataOpsReal<NDIM,double> ops;

    const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    if ( !pgeom->isTouchesBlockSingularBoundary() ) { // 不接触奇异点

        // 计算速度修正量
	ANCOnBox(patch,ghost_box,velocity_tmp);

        // 修正速度
        // v_scratch = v_current + v_tmp/4*anc
        double alpha = 0.25*d_anc;
#if ERROR_DIFFUSE
        ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_current, patch_box );
#else
        ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_scratch, patch_box );
#endif

    } else {
        // 获取奇异点信息
        hier::BoxArray<NDIM> fill_boxes;
        tbox::Array<int> location_indexes;
        tbox::Array<tbox::Array<tbox::Pointer<hier::Patch<NDIM > > > > singu_patches;
        tbox::Array< tbox::Array< int > > blocks;
	int codim = NDIM;
        const hier::Box<NDIM> &patch_box = patch.getBox();
        const hier::IntVector<NDIM> gcw(1);
	pgeom->getBlockInternalSingularBoundaryFillBoxes(
                                           fill_boxes, 
		                           location_indexes, 
					   singu_patches, 
					   blocks, 
					   codim, 
					   patch_box, 
					   gcw);

        // 处理非奇异点
        // (1) 计算速度修正量
        hier::BoxList<NDIM> box_list(ghost_box);
        box_list.removeIntersections(fill_boxes);
        hier::BoxArray<NDIM> boxes(box_list);
        for ( int i=0; i<boxes.size(); i++ ) {
  	    ANCOnBox(patch,boxes(i),velocity_tmp);
        }

        for ( int i=0; i<fill_boxes.size(); i++ ) {	
           velocity_tmp->fillAll(0.0, fill_boxes(i));  // 奇异点修正量清零
        }

        // (2) 修正速度
        double alpha = 0.25*d_anc;
#ifdef ERROR_DIFFUSE
        // v_scratch = v_current + v_tmp/4*anc
        ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_current, patch_box );
#else
        // v_scratch = v_scratch + v_tmp/4*anc
        ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_scratch, patch_box );
#endif


        // 处理奇异点
        for ( int i=0; i<fill_boxes.size(); i++ ) {	
           // (1) 计算速度修正量
  	   for(int n=0; n<singu_patches[i].size(); n++) {
	      tbox::Pointer<hier::Patch<NDIM> > singu_patch = singu_patches[i][n];
	      const hier::Box<NDIM> &singu_box = singu_patch->getBox();
              tbox::Pointer< pdat::NodeData<NDIM,double> > singu_velocity_tmp 
	         = new pdat::NodeData<NDIM,double>(singu_box,NDIM,1); //保存纠正量
              singu_velocity_tmp->fillAll(0.0);
	      ANCOnBox(*singu_patch,singu_patch->getBox(),singu_velocity_tmp);
              ops.add(velocity_tmp,velocity_tmp,singu_velocity_tmp,fill_boxes(i));
  	   }

           int n = singu_patches[i].size();
           double alpha = d_anc / n ;
#ifdef ERROR_DIFFUSE
           // v_scratch = v_current + v_tmp/n*anc
           ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_current, fill_boxes(i));
#else
           // v_scratch = v_scratch + v_tmp/n*anc
           ops.axpy(velocity_scratch, alpha, velocity_tmp, velocity_scratch, fill_boxes(i));
#endif
        }
    }
}

void Dynamics::ANCOnBox(hier::Patch<NDIM>& patch,
                        const hier::Box<NDIM> &box,
                        tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_tmp)
{
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_scratch =
      patch.getPatchData(d_velocity_scratch_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!velocity_scratch.isNull());
#endif 
   alternal_node_coupler_(
         ifirst_cells(0),ilast_cells(0),
         ifirst_cells(1),ilast_cells(1),
         ghost_cells(0),ghost_cells(1),
	 box.lower(0), box.upper(0),
	 box.lower(1), box.upper(1),
         velocity_scratch ->getPointer(),// 结点速度 
         velocity_tmp ->getPointer(),
	 d_xi);
}

/*
*************************************************************************
* 计算网格片内部网格结点的速度.				                *
*************************************************************************
*/
void Dynamics::computeVelocityOnPatch(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const double dt)
{
   (void) time;

   const hier::IntVector<NDIM> ghost_cells(1);

   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_scratch =
      patch.getPatchData(d_velocity_scratch_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity_new =
      patch.getPatchData(d_velocity_new_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > speedup =
      patch.getPatchData(d_speedup_scratch_id);

   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!velocity_scratch.isNull());
   assert(!velocity_new.isNull());
   assert(!speedup.isNull());
#endif 
   
   compute_velocity_(             
         ifirst_cells(0),ilast_cells(0),
         ifirst_cells(1),ilast_cells(1),
         ghost_cells(0),ghost_cells(1),
         velocity_scratch ->getPointer(),
         velocity_new     ->getPointer(),      // 结点速度 
         speedup  ->getPointer(),      // 结点加速度
         dt);

   postprocessVelocityOnPatch(patch,
                              time,
                              dt);
}

/*
*************************************************************************
* 后处理物理边界上的结点速度.                  	                       	*
*************************************************************************
*/

void Dynamics::postprocessVelocityOnPatch(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const double dt)
{
   (void) time;

   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates =
      patch.getPatchData(d_coordinates_scratch_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity_new_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coordinates.isNull());
   assert(!velocity.isNull());
#endif 

   // 是否接触物理边界?
   tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > patch_geometry =
                  patch.getPatchGeometry();

   const hier::Box<NDIM>& patch_box = patch.getBox();

   int btype=1;
   hier::BoxArray<NDIM> fill_boxes;
   tbox::Array<int> loc_indexes;
   tbox::Array<int> bdry_cond;
   patch_geometry->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, bdry_cond,
                                          btype, patch_box, ghost_cells);

   for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);

          // 处理边界条件
          if( bdry_cond[num] == SYMM_BDY2D_COND_TYPE ||
              bdry_cond[num] == WALL_BDY2D_COND_TYPE ||
              bdry_cond[num] == FREE_BDY2D_COND_TYPE ||   
              bdry_cond[num] == DIRI_BDY2D_COND_TYPE ||   
              bdry_cond[num] == OFLO_BDY2D_COND_TYPE ||   
              bdry_cond[num] == ZERO_BDY2D_COND_TYPE )  {   


                 postprocess_velocity_(
                    ifirst_cells(0),ilast_cells(0),
                    ifirst_cells(1),ilast_cells(1),
                    ghost_cells(0),ghost_cells(1),
                    fill_box.lower(0),fill_box.upper(0),
                    fill_box.lower(1),fill_box.upper(1),
                    coordinates->getPointer(),    // 结点坐标
                    velocity->getPointer(),      // 结点速度  
                    loc_indexes[num],bdry_cond[num]);
         }

   }
}


/*
*************************************************************************
* 计算网格片内部网格结点的坐标值、单元密度和能量.		        *
*************************************************************************
*/
void Dynamics::computeGridDensityEnergyOnPatch(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const double dt)
{
   (void) time;
   
   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates =
      patch.getPatchData(d_coordinates_current_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates_new =
      patch.getPatchData(d_coordinates_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > density =
      patch.getPatchData(d_density_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > density_new =
      patch.getPatchData(d_density_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > energy =
      patch.getPatchData(d_energy_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > energy_new =
      patch.getPatchData(d_energy_new_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > gamma  =
      patch.getPatchData(d_gamma_current_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
      patch.getPatchData(d_pressure_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > paq =
      patch.getPatchData(d_paq_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > viscosity =
      patch.getPatchData(d_viscosity_scratch_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity_new_id);

   const hier::IntVector<NDIM> ghosts(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();


   compute_griddensityenergy_(             
         ifirst_cells(0),ilast_cells(0),
         ifirst_cells(1),ilast_cells(1),
         ghosts(0),ghosts(1),
         coordinates     ->getPointer(),      
         coordinates_new ->getPointer(),     
         density         ->getPointer(),    
         density_new     ->getPointer(),   
         energy          ->getPointer(),  
         energy_new      ->getPointer(), 
         pressure        ->getPointer(),
         paq             ->getPointer(),
         velocity        ->getPointer(),
         viscosity       ->getPointer(),
         gamma           ->getPointer(),
         dt
	 );

}

/*
*************************************************************************
* 填充物理边界影像区.
*************************************************************************
*/
void Dynamics::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   if(intc_name!= "SIDE_PRESSURE" )return;

   tbox::Pointer< pdat::NodeData<NDIM,double> > coordinates =
      patch.getPatchData(d_coordinates_scratch_id);
   tbox::Pointer< pdat::NodeData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity_scratch_id);

   tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > paq =  
      patch.getPatchData(d_paq_scratch_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > mass =  
      patch.getPatchData(d_mass_scratch_id);

   const hier::IntVector<NDIM> ghost_cells(1);
   const hier::IntVector<NDIM> ifirst_cells = (patch.getBox()).lower();
   const hier::IntVector<NDIM> ilast_cells  = (patch.getBox()).upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coordinates.isNull());
   assert(!velocity.isNull());
   assert(!density.isNull());
   assert(paq->getGhostCellWidth()==ghost_cells);
   assert(mass->getGhostCellWidth()==ghost_cells);
#endif

   // 是否接触物理边界?
   tbox::Pointer< geom::BlockDeformingPatchGeometry<NDIM> > patch_geometry =
                  patch.getPatchGeometry();

   const hier::Box<NDIM>& patch_box = patch.getBox();

   for(int btype=1; btype<=NDIM; btype++) {

      hier::BoxArray<NDIM> fill_boxes;
      tbox::Array<int> loc_indexes;
      tbox::Array<int> bdry_cond;
      patch_geometry->getPhysicalBoundaryFillBoxes(fill_boxes,loc_indexes, bdry_cond,
                                          btype, patch_box, ghost_cells);

      for (int num=0; num < fill_boxes.size(); num++) {

         hier::Box<NDIM> fill_box = fill_boxes(num);

	 if(btype==1) {
            assert(bdry_cond[num] == SYMM_BDY2D_COND_TYPE ||
                   bdry_cond[num] == WALL_BDY2D_COND_TYPE ||
                   bdry_cond[num] == DIRI_BDY2D_COND_TYPE ||
                   bdry_cond[num] == ZERO_BDY2D_COND_TYPE ||
                   bdry_cond[num] == OFLO_BDY2D_COND_TYPE ||
                   bdry_cond[num] == FREE_BDY2D_COND_TYPE);
	 }
	 if(btype==2) {
            assert(bdry_cond[num] == X_SYMM_CORNER_TYPE ||
                   bdry_cond[num] == R_SYMM_CORNER_TYPE ||
                   bdry_cond[num] == X_WALL_CORNER_TYPE ||
                   bdry_cond[num] == R_WALL_CORNER_TYPE ||
                   bdry_cond[num] == X_FREE_CORNER_TYPE ||
                   bdry_cond[num] == R_FREE_CORNER_TYPE  
                   );
         }

	 physical_boundary_conditions_(
              ifirst_cells(0),ilast_cells(0),
              ifirst_cells(1),ilast_cells(1),
              ghost_cells(0),ghost_cells(1),
              fill_box.lower(0),fill_box.upper(0),
              fill_box.lower(1),fill_box.upper(1),
              coordinates->getPointer(),    // 结点坐标
              velocity->getPointer(),       // 结点速度
              density->getPointer(),        // 密度
              mass->getPointer(),        // 密度
              paq    ->getPointer(),        // 压力和人为粘力的总和.
              btype,                    // 边界box的类型
              loc_indexes[num],               // 边界box的位置
              bdry_cond[num]);                   // 边界条件
      }
   }
}

/*
*************************************************************************
* 将绘图量注册到JaVis数据输出器.                                        *
*************************************************************************
*/
void Dynamics::registerPlotData(
   tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(javis_writer.isNull()));
#endif

   if (!(javis_writer.isNull())) {

      javis_writer->registerNodeCoordinates(d_coordinates_current_id);
      
      javis_writer->registerPlotQuantity("Velocity",
                                           "VECTOR",
                                            d_velocity_current_id);
      javis_writer->registerPlotQuantity("Density",
                                           "SCALAR",
                                            d_density_current_id);
      javis_writer->registerPlotQuantity("Energy",
                                           "SCALAR",
                                            d_energy_current_id);

   } 

}

/*
*************************************************************************
* 将该对象的所有数据成员输出到指定的输出流.                             *
*************************************************************************
*/

void Dynamics::printClassData(ostream &os) const 
{

   os << "\nDynamics::printClassData..." << endl;
   os << "d_object_name = " << d_object_name << endl;

   os << "Parameters for physical problem ..." << endl;
   os << "   d_model = " << d_model << endl;

   os << "Numerical method description and ghost sizes..." << endl;
   os << "   d_anc = " << d_anc << endl;
   os << "   d_xi = " << d_xi << endl;
   os << "   d_a0 = " << d_a0 << endl;
   os << "   d_b0 = " << d_b0 << endl;

}


/*
*************************************************************************
* 从输入数据库读取数据,它们将覆盖从重启动数据库中读取的参数.            *
*************************************************************************
*/

void Dynamics::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if ( !is_from_restart ) {
      d_model   = db->getString("model");

      d_anc     = db->getDoubleWithDefault("anc",0.0);
      d_xi      = db->getDoubleWithDefault("xi",1.0);
      d_a0      = db->getDoubleWithDefault("a0",1.44);
      d_b0      = db->getDoubleWithDefault("b0",0.0);

      d_initial_gamma = db->getDouble("initial_gamma");
      d_initial_density  = db->getDouble("initial_density");
      d_initial_energy  = db->getDouble("initial_energy");
      d_initial_center_energy  = db->getDouble("initial_center_energy");
      db->getDoubleArray("initial_velocity",d_initial_velocity, NDIM);

   }
}


/*
*************************************************************************
* 从重启动数据库读入该对象的数据成员值 			*
*************************************************************************
*/
void Dynamics::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
   }

   d_anc     = db->getDouble("d_anc");
   d_xi      = db->getDouble("d_xi");
   d_a0      = db->getDouble("d_a0");
   d_b0      = db->getDouble("d_b0");
   d_model   = db->getString("d_model");

   d_initial_gamma = db->getDouble("d_initial_gamma");
}


/*
*************************************************************************
* 将该对象的数据成员输出到重启动数据库, 			*
*************************************************************************
*/
void Dynamics::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putString("d_model",d_model);
   db->putDouble("d_anc",d_anc);
   db->putDouble("d_xi",d_xi);
   db->putDouble("d_a0",d_a0);
   db->putDouble("d_b0",d_b0);

   db->putDouble("d_initial_gamma",d_initial_gamma);
}
