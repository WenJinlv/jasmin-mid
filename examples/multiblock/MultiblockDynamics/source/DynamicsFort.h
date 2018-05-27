//
// 文件  :	DynamicsFort.h
// 软件包:	JASMIN application
// 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:	$Revision: 1.14 $
// 修改  :	$Date: 2006/08/07 03:13:21 $
// 描述  :	求解流体力学方程的Fortran 77子程序的接口说明.
//

extern "C" {
  void initialize_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const double *, 
  double *,
  const double&);

  void zerobottom_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  double *);

  void zeroleft_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  double *);


  void zerobottom_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  double *);


  void stable_dt_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const double *, const double *, const double *, const double *,
  double &);

  void physical_boundary_conditions_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const int& , const int&, const int&, const int&, 
  const double *, const double *, const double *, const double *,
  double*,
  const int &, const int &, const int &);
 
  void preprocess_dynamics_state_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const double *, const double *, double *, 
  const double *, const double *, double *, double *, double *,
  const double *, const double&, const double&);

   void compute_side_pressure_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const double *, const double *,
  double *, double *);


  void compute_speedup_iga_(
  const int& , const int&, const int&, const int&, const int&, const int&, 
  const int&, const int&, const int&, const int&,
  const double *, const double *, 
  const double *, const double *, 
  double *, double *);
	     
  void postprocess_speedup_iga_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const int& , const int&, const int&,  const int&, const int&, const int&,
  const double *, const double *, const double *,
  double *, double *, const int &);
 
  void alternal_node_coupler_(
  const int& , const int&, const int&, const int&, const int&, const int&, 
  const int& , const int& ,const int& ,const int& ,
  double *, double *, const double &);

  void compute_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  double *, const double *, double *, const double & );

  void postprocess_velocity_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const int& , const int&, const int&,  const int&,
  const double *,  double *, 
  const int&, const int&);

  void compute_griddensityenergy_(
  const int& , const int&, const int&, const int&, const int&, const int&,
  const double * , double *,
  const double * , double *,
  const double * , 
  const double * , double *,
  const double * , const double * ,const double * ,
  const double * , const double& );

  void output_patch_data_dynamics_(
  const int&, const int&, const int&, const int&, const int&, const int&, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const int*, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *);

  void output_patch_data_dynamics0_(
  const int&, const int&, const int&, const int&, const int&, const int&, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *, 
  const int&, const int&, const double *);
}

