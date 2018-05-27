//
// 文件名: MD.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800 (浜, 28  9 2007) $
// 描述  : 分子动力学问题的网格片时间积分算法类的实现.
//

#ifndef included_MD_C
#define included_MD_C

#include "MD.h"

#include "tbox/PIO.h"
#include "tbox/MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"
#include "BoundaryBox.h"
#include "CellData.h"
#include "ParticleData.h"
#include "FaceData.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


#include "UniRectangularPatchGeometry.h"

// Fortran数值子程序接口.
#include "MDFort.h"

// 逻辑值常量.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

using namespace std;

/****************************************************************
* 构造函数.
*****************************************************************/
MD::MD(const string& object_name,
       tbox::Pointer<tbox::Database> input_db,
       tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!input_db.isNull());
    assert(!grid_geom.isNull());
#endif

    d_object_name = object_name;

    tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    d_grid_geometry = grid_geom;
    d_integrator_step = 0;
    d_particle_memory_grow = 1.2;

    // 缺省的影像单元宽度.
    d_nghosts = hier::IntVector<NDIM>(1);

    // 读取输入参数和重启动数据, 并用其初始化对象.
    bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart){
        getFromRestart();
    }
    getFromInput(input_db, is_from_restart);


    const double* dx  = d_grid_geometry->getDx();
    const double* xlo = d_grid_geometry->getXLower();
    const double* xhi = d_grid_geometry->getXUpper();
    hier::IntVector<NDIM> periodic_shift = d_grid_geometry->getPeriodicShift();
    const int myrank = tbox::MPI::getRank() ;


    hier::BoxArray<NDIM> aboxes = d_grid_geometry->getPhysicalDomain();
    if(!d_grid_geometry->getDomainIsSingleBox()){
        TBOX_ERROR(" physical domain is not a single box" << endl);
    }
    const hier::Index<NDIM> ifirst= aboxes.getBox(0).lower();
    const hier::Index<NDIM> ilast = aboxes.getBox(0).upper();


    /*
     * 将与问题相关的参数传递给Fortran程序.
     */
    setglobalparam_(&myrank,periodic_shift,dx,xlo,xhi,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM>2)
                    ifirst(2),ilast(2),
#endif
                    &d_left_obj_velocity,d_left_obj_position,
                    &d_right_obj_velocity,d_right_obj_position,
                    &d_triangle_vertex, &d_triangle_angle) ;


    // 创建所有变量及数据片索引号, 注册可视化数据片.
    registerModelVariables();

}

/****************************************************************
* 构造函数.
*****************************************************************/
MD::~MD()
{
    tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/
void MD::registerModelVariables()
{
#ifdef DEBUG_CHECK_ASSERTIONS
#endif

    hier::IntVector<NDIM> zero_ghosts(0);

    hier::VariableDatabase<NDIM>* variable_db =
        hier::VariableDatabase<NDIM>::getDatabase();
    d_current = variable_db->getContext("CURRENT");
    d_new = variable_db->getContext("NEW");

    d_plot_context = d_current;

    d_average_number_in_cell=1;
    d_Cu_particle_var   = new pdat::ParticleVariable<NDIM,double>
                       ("Cu_particle", NDIM*2, NDIM+2, d_average_number_in_cell);
    d_Cu_particle_id = variable_db->registerVariableAndContext(
       d_Cu_particle_var, d_new,  d_nghosts);

    d_particle_density_var  = new pdat::CellVariable<NDIM,double>("ParticleDensity", 1);
    d_particle_density_id = variable_db->registerVariableAndContext(
                                d_particle_density_var, d_current, zero_ghosts);
    d_patch_mapping_var  = new pdat::CellVariable<NDIM,double>("PatchMapping", 1);
    d_patch_mapping_id = variable_db->registerVariableAndContext(
                                d_patch_mapping_var, d_current, zero_ghosts);

/*
    d_number_particles_on_patch =
        tbox::Statistician::getStatistician()->
        getStatistic("number_particles_on_patch", "PATCH_STAT");
*/
}

double MD::getPatchDt(hier::Patch<NDIM>& patch,
                      const double  time,
                      const bool    initial_time,
                      const int     flag_last_dt,
                      const double  last_dt,
                      const string& intc_name)
{       
    return(d_time_step);
}


/*
*************************************************************************
*                                                                       *
* 注册JaVis数据输出器, 以采用JaVis工具对绘图文件进行后处理.           *
*                                                                       *
*************************************************************************
*/
void MD::registerJaVisDataWriter(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
#endif
    d_Javis_writer = viz_writer;

    if (!(d_Javis_writer.isNull())) {

        d_Javis_writer->registerParticleCoordinates( d_Cu_particle_id, 0);

        d_Javis_writer->registerPlotQuantity("ParticleDensity",
                                             "SCALAR", d_particle_density_id);

        d_Javis_writer->registerPlotQuantity("PatchMapping",
                                             "SCALAR", d_patch_mapping_id);
    }
}

/********************************************************************
*  初始化积分构件.
*********************************************************************/
void MD::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(intc!=NULL);
#endif

    const string intc_name = intc->getName();

    if(intc_name=="INIT") {  // 初值构件.
        intc->registerInitPatchData(d_particle_density_id);
        intc->registerInitPatchData(d_patch_mapping_id);
        intc->registerInitPatchData(d_Cu_particle_id);
    }else if(intc_name=="PARTICLE_COMM") { // 负责粒子通信的粒子量通信构件.
        intc->registerRefinePatchData(d_Cu_particle_id,
                                      d_Cu_particle_id);
    }else if(intc_name=="PARTICLE_COMP") { // 更新粒子位置的数值构件.
    }else if(intc_name=="PARTICLE_REDUCE") {  // 规约构件：计算粒子总数目.
    }else {
        TBOX_WARNING("\n::initializeComponent() : component "
                     << intc_name <<" is not matched. "<<endl);
    }
}


/********************************************************************
*  计算时间步长.
*********************************************************************/
void MD::computeOnPatch(hier::Patch<NDIM>& patch,
                        const double  time,
                        const double  dt,
                        const bool    initial_time,
                        const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   #endif


    if(intc_name=="PARTICLE_COMP") {
        updateParticleOnPatch(patch,time,dt);
    } else{
        TBOX_WARNING("\n::computePatch() : component is not matched. "<<endl);
    }
}

/********************************************************************
*  初始化数据片.
*********************************************************************/
void MD::initializePatchData(
    hier::Patch<NDIM>& patch,
    const double  time,
    const bool    initial_time,
    const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
    assert(intc_name=="INIT");
   #endif

    (void) time;
    if (initial_time) {
        const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> > pgeom =
            patch.getPatchGeometry();
        const double* dx  = pgeom->getDx();
        const double* xlo = pgeom->getXLower();
        const double* xhi = pgeom->getXUpper();

        tbox::Pointer< pdat::ParticleData<NDIM,double> > particle_curr  =
            patch.getPatchData(d_Cu_particle_id);

        tbox::Pointer< pdat::CellData<NDIM, double> > particle_density  =
            patch.getPatchData(d_particle_density_id );

#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!particle_curr.isNull());
        assert(!particle_density.isNull());
#endif
        hier::IntVector<NDIM> particle_gcells = particle_curr->getGhostCellWidth();
        hier::IntVector<NDIM> density_gcells = particle_density->getGhostCellWidth();

        const hier::Index<NDIM> ifirst=patch.getBox().lower();
        const hier::Index<NDIM> ilast =patch.getBox().upper();

        int real_num = 0;
        const int dbl_depth = particle_curr->getDblAttrDepth();
        const int int_depth = particle_curr->getIntAttrDepth();
        int max_num = particle_curr->getMaxNumberOfParticles();

        //估算网格片上的粒子个数
        particle_density->fillAll(0.0);
        int is_compute_density = 1;


        parinit_(
            dx,xlo,xhi,
            ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM>2)
            ifirst(2),ilast(2),
#endif
            particle_gcells(0),
            particle_gcells(1),
#if (NDIM>2)
            particle_gcells(2),
#endif
            &max_num, &real_num,
            particle_curr->getDblAttrPointer(), &dbl_depth,
            particle_curr->getIntAttrPointer(), &int_depth,
            particle_curr->getIndexParticlesPointer(),

            density_gcells(0),
            density_gcells(1),
#if (NDIM>2)
            density_gcells(2),
#endif
            particle_density->getPointer(), &is_compute_density
        );


        //统计网格片上的粒子个数
        double *ptr = particle_density->getPointer() ;
        int nsize = patch.getBox().size();
        double aa = nsize ;
        for(int i=0; i<nsize; i++)aa += ptr[i];
        //考虑影像区.
        double bb = 1.2*particle_curr->getGhostBox().size()
             /(particle_curr->getBox().size()+0.1) ;
        int num = int(  aa * bb );

        //调整粒子数据片的内存
        particle_curr->resizeAttributeArray( num , false); 
        max_num = particle_curr->getMaxNumberOfParticles();
        is_compute_density = 0;


        parinit_(
            dx,xlo,xhi,
            ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM>2)
            ifirst(2),ilast(2),
#endif
            particle_gcells(0),
            particle_gcells(1),
#if (NDIM>2)
            particle_gcells(2),
#endif
            &max_num, &real_num,
            particle_curr->getDblAttrPointer(), &dbl_depth,
            particle_curr->getIntAttrPointer(), &int_depth,
            particle_curr->getIndexParticlesPointer(),

            density_gcells(0),
            density_gcells(1),
#if (NDIM>2)
            density_gcells(2),
#endif
            particle_density->getPointer(), &is_compute_density
        );
        //按单元紧凑格式存储粒子
        particle_curr->arrangeNewParticlesByCell(real_num);

        // 计算单元粒子密度
        if(real_num<1)particle_density->fillAll(0.0);
        if(real_num>0)
        computedensity_(
                     ifirst(0),ilast(0),
                     ifirst(1),ilast(1),
#if (NDIM>2)
                     ifirst(2),ilast(2),
#endif
                     particle_gcells(0),
                     particle_gcells(1),
#if (NDIM>2)
                     particle_gcells(2),
#endif
                     particle_curr->getIndexParticlesPointer(),
                     density_gcells(0),
                     density_gcells(1),
#if (NDIM>2)
                     density_gcells(2),
#endif
                     particle_density->getPointer()
                    );

      // record the number of particles in each patch
      /*
      int patch_number9 = patch.getPatchNumber();
      d_number_particles_on_patch->recordPatchStat(patch_number9,
        real_num * 1.0, 0);  
      */
   
      // record processor mapping of each patch
      const int myrank = tbox::MPI::getRank() ;
      tbox::Pointer< pdat::CellData<NDIM, double> > patch_mapping  =
        patch.getPatchData(d_patch_mapping_id );
       patch_mapping->fillAll(myrank*1.0);
 }

}

/********************************************************************
*  设置积分时间步数.
*********************************************************************/

void MD::setIntegratorStep(
    const int  integrator_step)
{
    d_integrator_step = integrator_step;
}
/********************************************************************
*  计算网格片中单元上的负载.
*********************************************************************/
void MD::computeWorkloadOnPatch(
    hier::Patch<NDIM>& patch,
    const double measured_load,
    const int    numeric_step,
    const int    workload_index)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   #endif
    const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> >
    pgeom = patch.getPatchGeometry();
    const double* dx  = pgeom->getDx();
    const double* xlo = pgeom->getXLower();
    const double* xhi = pgeom->getXUpper();
    const hier::Index<NDIM> ifirst=patch.getBox().lower();
    const hier::Index<NDIM> ilast =patch.getBox().upper();

    tbox::Pointer< pdat::CellData<NDIM,double> > load_weight  =
        patch.getPatchData(workload_index);
    hier::IntVector<NDIM> load_gcells = load_weight->getGhostCellWidth();
    hier::IntVector<NDIM> particle_gcells(0);

    load_weight->fillAll(0.0);

    // initial time

    if(d_integrator_step==0){
        const int real_num = 0;
        const int dbl_depth = 1;
        const int int_depth = 1;
        const int max_num = 1;
        tbox::Array<double> dbl_attr(10);
        tbox::Array<int>    int_attr(10);
        tbox::Array<int>    index_p(10);

        int is_compute_density = 1;
        parinit_(
            dx,xlo,xhi,
            ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM>2)
            ifirst(2),ilast(2),
#endif
            particle_gcells(0),
            particle_gcells(1),
#if (NDIM>2)
            particle_gcells(2),
#endif
            &max_num, &real_num,
            dbl_attr.getPointer(), &dbl_depth,
            int_attr.getPointer(), &int_depth,
            index_p.getPointer(),

            load_gcells(0),
            load_gcells(1),
#if (NDIM>2)
            load_gcells(2),
#endif
            load_weight->getPointer(), &is_compute_density
        );
    }
    else
    {

        tbox::Pointer< pdat::ParticleData<NDIM,double> > particle_new  =
            patch.getPatchData(d_Cu_particle_id);
        hier::IntVector<NDIM> ghost_cells = particle_new->getGhostCellWidth();

        double ttt = measured_load ;
        computeload_(&ttt, dx, xlo, xhi,
                     ifirst(0),ilast(0),
                     ifirst(1),ilast(1),
#if (NDIM>2)
                     ifirst(2),ilast(2),
#endif
                     ghost_cells(0),
                     ghost_cells(1),
#if (NDIM>2)
                     ghost_cells(2),
#endif
                     load_gcells(0),
                     load_gcells(1),
#if (NDIM>2)
                     load_gcells(2),
#endif
                     particle_new->getIndexParticlesPointer(),
                     load_weight->getPointer()
                    );
    }
}

/********************************************************************
*  计算粒子受力，更新粒子位置.
*********************************************************************/
void MD::updateParticleOnPatch
(hier::Patch<NDIM>& patch,
 const double current_time,
 const double dt)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   #endif
    const tbox::Pointer<geom::UniRectangularPatchGeometry<NDIM> >
    pgeom = patch.getPatchGeometry();
    const double* dx  = pgeom->getDx();
    const double* xlo = pgeom->getXLower();
    const double* xhi = pgeom->getXUpper();

    tbox::Pointer< pdat::ParticleData<NDIM,double> > particle_new  =
        patch.getPatchData(d_Cu_particle_id);
    double time8 = dt ;
    const hier::Index<NDIM> ifirst=patch.getBox().lower();
    const hier::Index<NDIM> ilast =patch.getBox().upper();
    hier::IntVector<NDIM> ghost_cells = particle_new->getGhostCellWidth();

    const int dbl_depth = particle_new->getDblAttrDepth();
    const int int_depth = particle_new->getIntAttrDepth();

    int max_num  = particle_new->getMaxNumberOfParticles();
    int real_num = particle_new->getNumberOfParticles();

      const int myrank = tbox::MPI::getRank() ;
      tbox::Pointer< pdat::CellData<NDIM, double> > patch_mapping  =
        patch.getPatchData(d_patch_mapping_id );
       patch_mapping->fillAll(myrank*1.0);
    if(real_num>0) {
        tbox::Array<double> auxi_pf(real_num*NDIM);
        updateparticles_(
            &time8,dx,xlo,xhi,
            ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM>2)
            ifirst(2),ilast(2),
#endif
            ghost_cells(0),
            ghost_cells(1),
#if (NDIM>2)
            ghost_cells(2),
#endif
            &max_num, &real_num,
            particle_new->getDblAttrPointer(), &dbl_depth,
            particle_new->getIntAttrPointer(), &int_depth,
            particle_new->getIndexParticlesPointer(),
            auxi_pf.getPointer()
        );
    }

// compute the number of particles in each cell
    tbox::Pointer< pdat::CellData<NDIM, double> > particle_density  =
          patch.getPatchData(d_particle_density_id );
    hier::IntVector<NDIM> density_gcells = particle_density->getGhostCellWidth();
//    if(real_num<1)particle_density->fillAll(0.0);
    particle_density->fillAll(0.0);
    if(real_num>0)
        computedensity_(
                     ifirst(0),ilast(0),
                     ifirst(1),ilast(1),
#if (NDIM>2)
                     ifirst(2),ilast(2),
#endif
                     ghost_cells(0),
                     ghost_cells(1),
#if (NDIM>2)
                     ghost_cells(2),
#endif
                     particle_new->getIndexParticlesPointer(),
                     density_gcells(0),
                     density_gcells(1),
#if (NDIM>2)
                     density_gcells(2),
#endif
                     particle_density->getPointer()
                    );

}

/********************************************************************
*  填充物理边界条件.
*********************************************************************/
void MD::setPhysicalBoundaryConditions(
    hier::Patch<NDIM>& patch,
    const double fill_time,
    const hier::IntVector<NDIM>& ghost_width_to_fill,
    const string& intc_name)
{
}

/********************************************************************
*  输出数据成员到重启动数据库.
*********************************************************************/
void MD::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

}

/********************************************************************
*  打印数据成员.
*********************************************************************/
void MD::printClassData(ostream& os) const
{


    os << "\nMD::printClassData..." << endl;
    os << "MD: this = " << (MD*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_grid_geometry = "
    << (geom::UniRectangularGridGeometry<NDIM>*)d_grid_geometry << endl;

    os<< "list (variable, context, patch data index) :\n"
    << "   " << d_Cu_particle_var->getName() <<"  " << d_new->getName()
    <<"  " << d_Cu_particle_id << "\n"
    << endl;
}

/********************************************************************
*  从输入文件读取数据成员. 
*********************************************************************/
void MD::getFromInput(tbox::Pointer<tbox::Database> db,
                      bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    d_time_step = 0.001;
    if (db->keyExists("time_step_length")) {
        d_time_step = db->getDouble("time_step_length");
    }
    if (db->keyExists("ghost_cells")) {
        int* tmp_nghosts = d_nghosts;
        db->getIntegerArray("ghost_cells", tmp_nghosts, NDIM);
    }


    //设置缺省值: 速度和相对位置.
    d_left_obj_velocity = 10.0;
    d_left_obj_position[0] = 0.2 ;
    d_left_obj_position[1] = 0.4 ;
    //设置缺省值: 速度和相对位置.
    d_right_obj_velocity = 10.0;
    d_right_obj_position[0] = 0.5 ;
    d_right_obj_position[1] = 0.7 ;
    //从输入文件中获取速度和相对位置.
    if (db->keyExists("left_obj_velocity"))
        d_left_obj_velocity = db->getDouble("left_obj_velocity");
    if (db->keyExists("left_obj_position"))
        db->getDoubleArray("left_obj_position",d_left_obj_position, 2);
    if (db->keyExists("right_obj_velocity"))
        d_right_obj_velocity = db->getDouble("right_obj_velocity");
    if (db->keyExists("right_obj_position"))
        db->getDoubleArray("right_obj_position",d_right_obj_position, 2);

    if (db->keyExists("triangle_vertex"))
        d_triangle_vertex = db->getDouble("triangle_vertex");
    if (db->keyExists("triangle_angle"))
        d_triangle_angle = db->getDouble("triangle_angle");

}

/********************************************************************
*  从输入文件读取数据成员. 
*********************************************************************/
void MD::getFromRestart()
{
#ifdef DEBUG_CHECK_ASSERTIONS
#endif
}



/********************************************************************
 *  规约构件：计算粒子总数目.
*********************************************************************/
void MD::reduceOnPatch( double *vector,
                             int     len,
                             hier::Patch<NDIM> &patch,
                             const double time,
                             const double dt,
                             const string& intc_name )
{
    tbox::Pointer< pdat::ParticleData<NDIM,double> > particle_new  =
        patch.getPatchData(d_Cu_particle_id);
    vector[0] += particle_new->getNumberOfParticles( particle_new->getBox() );
}

#endif
