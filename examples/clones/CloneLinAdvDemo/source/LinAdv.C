//
// 文件名: LinAdv.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 $
// 描述  : 线性对流问题的网格片时间积分算法类的实现.
//

#include "LinAdv.h"
#include "LinAdvFort.h"

#include "CellData.h"
#include "SideData.h"
#include "UniRectangularPatchGeometry.h"
#include "CommunicatorDatabase.h"

using namespace JASMIN;

/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
LinAdv::LinAdv(const string& object_name,
               tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = object_name;

   // 从输入文件读入X方向的对流速度.
   d_x_velocity = input_db->getDouble("constant_x_velocity");
   assert( d_x_velocity > 0.0 );  

   // 注册变量和数据片.
   registerModelVariables();
}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
LinAdv::~LinAdv()
{
};


/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void LinAdv::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   tbox::Pointer<pdat::CellVariable<NDIM,double> > d_uval 
       = new pdat::CellVariable<NDIM,double>("uval",1);  // 中心量，深度为1.
   tbox::Pointer<pdat::SideVariable<NDIM,double> > d_flux  
       = new pdat::SideVariable<NDIM,double>("flux",1);  // 边心量，深度为1.
  
   // 当前值上下文, 新值上下文, 演算上下文.
   tbox::Pointer<hier::VariableContext> d_current 
       = variable_db->getContext("CURRENT");
   tbox::Pointer<hier::VariableContext> d_new
       = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   tbox::Pointer<hier::VariableContext> d_record
       = hier::VariableDatabase<NDIM>::getDatabase()->getContext("RECORD");

   // 数据片<uval,current>: 存储当前时刻的 u 值, 影像区宽度为1.
   d_uval_current_id = variable_db->registerVariableAndContext(d_uval, d_current, 1);
   d_uval_record_id = variable_db->registerVariableAndContext(d_uval, d_record, 1);

   // 数据片<uval,new>: 存储新时刻的 u 值, 影像区宽度为0.
   d_uval_new_id = variable_db->registerVariableAndContext(d_uval, d_new);

   // 数据片<flux,new>: 存储新时刻的 f 值, 影像区宽度为0.
   d_flux_new_id = variable_db->registerVariableAndContext(d_flux, d_new);

}


/*************************************************************************
 *                                                                        
 * 注册绘图量.
 *                                                                  
 *************************************************************************/
void LinAdv::registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
   javis_writer-> registerPlotQuantity("U", "SCALAR", d_uval_current_id);
   javis_writer-> registerPlotQuantity("U_REC", "SCALAR", d_uval_record_id);
}


/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void LinAdv::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
   const string &intc_name = intc->getName();

   if ( intc_name=="ALLOC_NEW_DATA") {  // 内存构件.
         intc->registerPatchData(d_uval_new_id);
         intc->registerPatchData(d_flux_new_id);

    } else if ( intc_name=="INIT_SET_VALUE") {  // 初值构件. 
         intc->registerInitPatchData(d_uval_current_id);
         intc->registerInitPatchData(d_uval_record_id);

    } else if(  intc_name=="STEP_SIZE" ) {  // 步长构件.

    } else if(  intc_name=="INIT_POST_SOLUTION" ) {  // 数值构件.

    } else if(  intc_name=="POST_SOLUTION" ) {  // 数值构件.

    } else if(  intc_name=="PREP_SOLUTION" ) {  // 数值构件.

    } else if ( intc_name=="COMPUTE") {  // 数值构件: 采用迎风格式，需要u的影像区数据.
          intc->registerRefinePatchData(d_uval_current_id, d_uval_current_id);

    } else if ( intc_name=="COPY_SOLUTION") {  // 复制构件.
          intc->registerCopyPatchData(d_uval_current_id, d_uval_new_id);   

    } else if ( intc_name=="BCAST_SOLUTION") {  // 广播构件.
          intc->registerPatchData(d_uval_current_id);   

    } else if ( intc_name=="INIT_COLL_SOLUTION") {  // 初始时刻的汇集构件.
          intc->registerPatchData(d_uval_current_id);   

    } else if ( intc_name=="REDUCTION_SOLUTION") {  // 汇集构件.
          intc->registerPatchData(d_uval_new_id);   

    } else if ( intc_name=="PULSATION_SOLUTION") {  // 脉动构件.
          intc->registerPatchData(d_uval_current_id);   

    } else if ( intc_name=="INIT_FILL_GHOST") {  // 汇集构件.
          intc->registerRefinePatchData(d_uval_current_id, d_uval_current_id);

    } else {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
    }

}

/*************************************************************************
 *
 *  初始化数据片（支持初值构件）.
 *
 ************************************************************************/
void LinAdv::initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name)
{
   if (initial_time) {
      tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
                                    patch.getPatchData(d_uval_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > uval_record =
                                    patch.getPatchData(d_uval_record_id);

      const hier::IntVector<NDIM> &ghost_cells = uval_current->getGhostCellWidth();

      const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const hier::Index<NDIM> &ifirst=patch.getBox().lower();
      const hier::Index<NDIM> &ilast =patch.getBox().upper();

      const int clone_number = pgeom->getCloneNumber();
      initsetvalue_(dx,xlo,clone_number,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    uval_current->getPointer(),
                    uval_record->getPointer());
   }

}

/*************************************************************************
 *
 *  计算稳定性时间步长（支持步长构件）.
 *
 ************************************************************************/
double LinAdv::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom
                       = patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();

   return 0.8*dx[0]/d_x_velocity;

}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 *
 * 该函数基于显式迎风格式，实现通量和守恒量的计算.
 ************************************************************************/
void LinAdv::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{
   // 获取几何信息
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* xlo = pgeom->getXLower();

   // 获取网格片索引范围
   const hier::Index<NDIM> &ifirst=patch.getBox().lower();
   const hier::Index<NDIM> &ilast =patch.getBox().upper();

   // 获取数据片
   tbox::Pointer< pdat::CellData<NDIM,double> > uval_current=
                                    patch.getPatchData(d_uval_current_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > uval_record =
                                    patch.getPatchData(d_uval_record_id);

   hier::IntVector<NDIM> ghost_cells = uval_current->getGhostCellWidth();

   // 克隆号.
   const int number_clones = pgeom->getNumberClones();
   const int clone_number = pgeom->getCloneNumber();

   if(intc_name=="COMPUTE") {

      tbox::Pointer< pdat::CellData<NDIM,double> > uval_new =
                                    patch.getPatchData(d_uval_new_id);
      tbox::Pointer< pdat::SideData<NDIM,double> > flux_new =
                                    patch.getPatchData(d_flux_new_id);

      // 计算通量
      advancepatch_(d_x_velocity, dt, dx,xlo,
                    ifirst(0),ilast(0),
	    	    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
	 	    ghost_cells(0),
		    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
		    uval_current->getPointer(),
		    flux_new->getPointer(0),
		    flux_new->getPointer(1)
#if (NDIM==3)
                    ,flux_new->getPointer(2)
#endif
                    );

      // 更新守恒量
      conspatch_(dx, xlo,
                 ifirst(0),ilast(0),
                 ifirst(1),ilast(1),
#if (NDIM==3)
                 ifirst(2),ilast(2),
#endif
                 ghost_cells(0),
                 ghost_cells(1),
#if (NDIM==3)
                 ghost_cells(2),
#endif
	         uval_current->getPointer(),
	         flux_new->getPointer(0),
	         flux_new->getPointer(1),
#if (NDIM==3)
                 flux_new->getPointer(2),
#endif
	         uval_new->getPointer() );

   } else if(intc_name=="INIT_POST_SOLUTION") {

      initpostsolution_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    number_clones,
	            uval_current->getPointer(),
	            uval_record->getPointer());

   } else if(intc_name=="POST_SOLUTION") {

      postsolution_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    number_clones,
	            uval_current->getPointer(),
	            uval_record->getPointer());

   } else if(intc_name=="PREP_SOLUTION") {

      prepsolution_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
                    clone_number,
	            uval_current->getPointer());

   } else if(intc_name=="INIT_FILL_GHOST") {};

}

/*************************************************************************
 * 
 ************************************************************************/
void LinAdv::pulseOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const int     src_clone_number,
                          const string& intc_name)
{
   // 获取几何信息
   const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom = 
                                    patch.getPatchGeometry();

   // 获取网格片索引范围
   const hier::Index<NDIM> &ifirst=patch.getBox().lower();
   const hier::Index<NDIM> &ilast =patch.getBox().upper();

   // 克隆号.
   const int clone_number = pgeom->getCloneNumber();

   // 仅在属主网格层上执行.
   if(intc_name=="PULSATION_SOLUTION") {

      tbox::Pointer< pdat::CellData<NDIM,double> > uval_current =
                                    patch.getPatchData(d_uval_current_id);
      tbox::Pointer< pdat::CellData<NDIM,double> > uval_record =
                                    patch.getPatchData(d_uval_record_id);
      hier::IntVector<NDIM> ghost_cells = uval_current->getGhostCellWidth();

      // 元素值*(src_clone_number+1)
      pulseclonecomputing_(
                    clone_number,
                    src_clone_number,
                    ifirst(0),ilast(0),
	    	    ifirst(1),ilast(1),
#if (NDIM==3)
                    ifirst(2),ilast(2),
#endif
	 	    ghost_cells(0),
		    ghost_cells(1),
#if (NDIM==3)
                    ghost_cells(2),
#endif
		    uval_current->getPointer(),
		    uval_record->getPointer()
                    );

   }

}

/*************************************************************************
 *
 *  填充物理边界条件.
 *
 ************************************************************************/
void LinAdv::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   // 获取数据片及其影像区宽度
   tbox::Pointer< pdat::CellData<NDIM,double> > uval_current
                          = patch.getPatchData(d_uval_current_id);
   const hier::IntVector<NDIM> &ghost_cells = uval_current->getGhostCellWidth();

   // 获取网格片索引范围
   const hier::Box<NDIM> &patch_box = patch.getBox();
   const hier::Index<NDIM> &ifirst=patch.getBox().lower();
   const hier::Index<NDIM> &ilast =patch.getBox().upper();

   // 获取(棱型)物理边界影像区
   const tbox::Pointer< hier::PatchGeometry<NDIM> > pgeom =
                                       patch.getPatchGeometry();
   hier::BoxArray<NDIM> fill_boxes;
   tbox::Array<int> loc_indexes;
   int btype=1;
   pgeom->getPhysicalBoundaryFillBoxes(fill_boxes,
                                  loc_indexes,
                                  btype, 
                                  patch_box,
                                  ghost_cells);

   const int clone_number  = pgeom->getCloneNumber();

   // 遍历影像区
   for (int i=0; i < fill_boxes.size(); i++) {
        // 获取物理边界影像区索引范围
        const hier::Index<NDIM> &igfirst=fill_boxes(i).lower();
	const hier::Index<NDIM> &iglast =fill_boxes(i).upper();

        // 调用FORTRAN函数，根据物理边界条件，填充影像区数据
        fillphysicbdry_(
                ifirst(0),ilast(0),
                ifirst(1),ilast(1),
#if (NDIM==3)
                ifirst(2),ilast(2),
#endif
                igfirst(0),iglast(0),
                igfirst(1),iglast(1),
#if (NDIM==3)
                igfirst(2),iglast(2),
#endif
                loc_indexes[i],
		ghost_cells(0),
		ghost_cells(1),
#if (NDIM==3)
                ghost_cells(2),
#endif
                clone_number,
		uval_current->getPointer() );
   }
}


