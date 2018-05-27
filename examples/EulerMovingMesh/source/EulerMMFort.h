//
// 文件名: EulerMMFort.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.2 $
// 修改  : $Date: 2007/01/06 09:07:25 $
// 描述  : F77外部子程序接口.
//

#include <math.h>
#include <signal.h>

extern "C" {
 
  void setinitialvalue_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const int& , const double&, double*);
 
  void stabledt_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double& , 
  const double*,
  const double*,
  const double*,
  double&);
 
  void solvepde_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double& , const double& , 
  const double*,
  const double*,
  const double*,
  const double*,
  const double*,
  double*,
  double*,
  double*);
 
  void computegradient_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double*,
  const double*,
  const double*,
  double*,
  double*);
 
  void computewcnew_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double&,const double&,
  const double&,const double&,
  const double*,
  const double*,
  const double*,
  double*);
 
  void wcnewiteration_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double*,
  double*);

  void movemesh_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const double*,
  const double*,
  const double*,
  double*, double*);

  void movesolution_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const double*,
  const double*,
  const double*,
  const double*,
  const double*,
  double*,
  const double*,
  const double*,
  double*,
  double*);

  void setphysbdryforcells_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& , const int& ,
  double*);

  void setphysbdryfornodes_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const double& , const double&, const double&,
  const double& , const double&, const double&,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  double*, double*);

  void postprocessbdrymesh_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const double*, const double*,
  double *, double *);

}
