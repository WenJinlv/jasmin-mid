//
// 文件名: LinAdvFort.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.2 $
// 修改  : $Date: 2007/01/06 09:07:25 $
// 描述  : F77外部子程序接口.
//

#include <math.h>
#include <signal.h>

extern "C" {
 
  void initsetvalue_(
  const double *, const double *,
  const int& , const int& , const int& , const int& , 
#if (NDIM>2) 
  const int& , const int& ,
#endif 
  const int& , const int& , 
#if (NDIM>2) 
  const int& ,
#endif
  double*);

  void advancepatch_(
  const double&, const double&, const double *, const double *,
  const int& , const int& , const int&, const int& ,
#if (NDIM>2) 
  const int& , const int& ,
#endif 
  const int& , const int& ,
#if (NDIM>2) 
  const int& ,
#endif
  const double*, double *, double *
#if (NDIM>2)
  , double *
#endif
  );

  void conspatch_(
  const double *, const double *,
  const int& , const int& , const int& , const int& , 
#if (NDIM>2) 
  const int& , const int& ,
#endif 
  const int& , const int& ,
#if (NDIM>2) 
  const int& ,
#endif
  const double*, const double *, const double *, 
#if (NDIM>2)
  const double *,
#endif
  double *);

}
