//
// 文件名: SnSweep.h
// 软件包: JASMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 250 $
// 修改  : $Date: 2007-09-28 13:11:58 $  
// 描述  : Sn输运计算的网格片策略类实现.
//

#ifndef included_SnSweep
#define included_SnSweep

#include "JASMIN_config.h"

#include "StandardComponentPatchStrategy.h"
#include "SweepingPatchStrategy.h"
#include "UniRectangularGridGeometry.h"
#include "JaVisDataWriter.h"

using namespace JASMIN;

/**
 * @brief 该类基于标准构件网格片策略类 algs::StandardComponentPatchStrategy 
 * 和扫描计算网格片策略类 algs::SweepingPatchStrategy, 
 * 实现Sn中子输运方程的网格片计算子程序.
 * 
 * 另外, 该类使用了群结点量数据片, 来存储中子的多群角通量.
 *
 * 该类对象从输入文件读取如下参数:
 *    - model_name \n
 *      字符串, 模型名称.
 *    - matter \n
 *      字符串, 材料代号.
 *    - ngroup \n
 *      整型, 能群总数.
 *    - v  \n
 *      双精度浮点型数组, 能群系数.
 *    - MatterXXX \n
 *      数据库型, XXX为材料代码.
 *      - SGMt \n
 *        双精度浮点型数组, 输运截面.
 *      - SGMs \n
 *        双精度浮点型数组, 散射截面.
 *    - s_n \n
 *      整型, Sn数.
 *    - SYYY \n
 *      数据库型, YYY为Sn数.
 *      - u \n
 *        双精度浮点型数组, 方向余弦离散值.
 *      - w \n
 *        双精度浮点型数组, 方向权重.
 */
class SnSweep :
   public algs::StandardComponentPatchStrategy<NDIM>,
   public algs::SweepingPatchStrategy<NDIM>,
   public appu::JaVisDerivedDataStrategy<NDIM>
{
public:

   /*! @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示对象名称.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    * @param grid_geom   输入参数, 指针, 指向均匀矩形网格几何.
    */
   SnSweep(const string& object_name,
          tbox::Pointer<tbox::Database> input_db,
          tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geom);

   /*! @brief 析构函数.
    */
   virtual ~SnSweep();

   ///@name 重载类 algs::StandardComponentPatchStrategy 的成员函数
   //@{
   /*!
    * @brief 初始化指定的积分构件, 
    * @param intc 指针, 输入参数, 指向待初始化的积分构件对象.
    */
    void initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const;

   /**
    * @brief 当前积分构件为初值构件. 该函数初始化数据片.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param time  输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param intc_name 输入参数, 字符串, 表示当前积分构件的名称.
    *
    */
   void initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name);

   /*! 
    * @brief 当前积分构件为步长构件. 该函数计算时间步长.
    * @note 
    * 该程序没有使用步长构件, 所以该函数不会被调用.
    */
   double getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name){return -1.;};

   /*!
    * @brief 当前积分构件为数值构件. 该函数完成数值计算.
    */
   void computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name);

   /*!
    * @brief 当前积分构件为规约构件. 该函数计算单个网格片上的迭代误差.
    */
   void reduceOnPatch(
                    double *vector,
                    int     len,
                    hier::Patch<NDIM> &patch,
                    const double time,
                    const double dt,
                    const string& intc_name );

   /**
    * @brief 在积分构件中, 为数据片填充物理边界条件.
    *
    * @param patch     输入参数, 网格片类, 表示网格片.
    * @param fill_time 输入参数, 双精度浮点型, 表示当前时刻.
    * @param ghost_width_to_fill 输入参数, 整型向量, 表示物理边界影像区被填充的宽度.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    */
   void setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);
   //@}

   ///@name 重载类 algs::SweepingPatchStrategy 的成员函数
   //@{

   /*!
    *@brief 返回扫描方向个数.
    */ 
   int getNumberOfAngles( const string& intc_name ){return d_ms;};

   /*!
    *@brief 为扫描构件计算单元间的依赖关系.
    */ 
   void setCellsDataDependencyOnPatch(
                JASMIN::pdat::SideData<NDIM, int>& cell_depend,
                hier::Patch<NDIM> &patch,
		const double time,
		const double dt,
		const string& intc_name );

   /*!
    * @brief 为扫描构件完成指定单元序列的扫描计算.
    *
    * @param patch 输入参数, 网格片类, 表示网格片.
    * @param angle_id 输入参数, 整型，方向编号.
    * @param ncells   输入参数, 整型，待扫描单元的数目.
    * @param cells    输入参数, 整型数组，待扫描的单元索引.
    * @param time  输入参数, 双精度浮点型, 表示当前时刻.
    * @param dt    输入参数, 双精度浮点型, 表示时间步长.
    * @param cell_dependency    输入参数, 边心数据片, 存储单元间依赖关系.
    * @param intc_name 输入参数, 字符串, 表示积分构件的名称.
    */
   void sweepingOnPatch(hier::Patch<NDIM>& patch,
                        const int angle_id,
			const int ncells,
			const int cells[][NDIM],
			const double time, 
			const double  dt,
			const JASMIN::pdat::SideData<NDIM, int>& cell_depend,
			const string& intc_name );

   //@}

   ///@name 可视化相关函数
   //@{

   /*!
    * @brief 注册绘图量到JaVis数据输出器.
    * @param javis_writer 输入参数, 指针, 指向 JaVis 数据输出器.
    */
   void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer);

   /**
    * @brief 由基本量计算导出量并将其打包到缓冲区.
    * @param buffer 	输入参数, 指针, 指向双精度浮点型缓冲区.
    * @param patch 	输入参数, 网格片类, 网格片对象.
    * @param region 	输入参数, Box类, 待打包数据的索引区域.
    * @param variable_name 输入参数, 字符串, 导出量的名称.
    * @param depth_id   输入参数, 整型, 出量分量索引号.
    * @return  逻辑型, 真值表明导出量在该网格片上已被准确地定义.
    */
   bool packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch<NDIM>& patch,
      const hier::Box<NDIM>& region,
      const string& variable_name,
      int depth_id);
   //@}

   ///@name 其他
   //@{

#ifdef TESTING 
   /* @brief 设置当前时间步数 */
   void setTimeStepNumber(const int nsteps);

   /* @brief 设置精确解系数yd0,yd1 */
   void setSolutionCoeff(const double dt);

   void diffResultOnPatch(hier::Patch<NDIM> &patch);
#endif

   //@}

private:

   /*!
    * @brief 从输入数据库和重启动数据库读入数据.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);
  
   /*!
    * @brief 创建变量和数据片索引号.
    */
   void registerModelVariables();

   void computeScatterItemOnPatch(hier::Patch<NDIM>& patch);
   void computeSoureceItemOnPatch(hier::Patch<NDIM>& patch,const double dt);
   void updateVolumeFluxOnPatch(hier::Patch<NDIM>& patch);


   string d_object_name;         /*!< 对象名 */
   tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > d_grid_geometry; /*!< 网格几何 */
   tbox::Array<double> d_angles;
   
   ///@name 模型参数
   //@{
   string d_matter_code;
   int    d_mg;
   int    d_sn_n;
   int    d_ms;
   tbox::Array<double> d_vn;
   tbox::Array<double> d_sgm_t0;
   tbox::Array<double> d_sgm_s0;
   tbox::Array<double> d_sn_cos;
   tbox::Array<double> d_sn_weight;
   //@}


   ///@name 数据片索引号. 
   //@{
   int d_fi_cur_id;
   int d_fi_new_id;
   int d_wf_cur_id;
   int d_wf_new_id;
   int d_wf_scr_id;
   int d_sgm_t_id;
   int d_sgm_s_id;
   int d_qs_scr_id;
   int d_qnw_scr_id;
   //@}

   ///@name 其他
   //@{
   double d_xleng[NDIM];
   
#ifdef TESTING   
   double d_yd0, d_yd1;
   tbox::Array<double> d_ysgm_t, d_ysgm_s;
   int d_nsteps;
#endif

   //@}

};

#endif
