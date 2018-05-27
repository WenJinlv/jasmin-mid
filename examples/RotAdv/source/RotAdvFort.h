//
// 文件名: RotAdvFort.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.2 $
// 修改  : $Date: 2007/01/06 09:07:25 $
// 描述  : F77外部子程序接口.
//

#include <math.h>
#include <signal.h>

extern "C" {
 
  void initrotation_(
  const int& ,
  const double*, const double*,
  const double*, const double&,
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int&,
#if (NDIM>1)
  const int& ,  
#endif
#if (NDIM>2)
  const int& ,  
#endif
  double*);
 
  void stabledt_(
#if (NDIM==3)
  const double*, 
#endif
  const double*, const double*,
  const int& , const int& ,
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& ,
#if (NDIM>1)
  const int& ,
#endif
#if (NDIM>2)
  const int& ,
#endif
  const double*, 
  double&);
 
  void inittraceflux_(
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*, 
#if (NDIM>1)
  double*      , double*      , double*      , 
#endif
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      ); 
 
  void chartracing0_(
  const double&, const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& , const double*, const double*, 
#if (NDIM>2) 
  const double*, 
#endif 
  const int& ,
  const double*,
  double*      , double*      ,
  double*      , double*      , 
  double*      , double*);
 
#if (NDIM>1)
  void chartracing1_(
  const double&, const int& , const int& , const int& ,const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int& , const double*, const double*, 
#if (NDIM>2) 
  const double*, 
#endif 
  const int& , 
  const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*);
 
#if (NDIM>2)
  void chartracing2_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const int& , const double*, const double*, const double *, const int& ,
  const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*);
#endif
#endif
 
  void fluxcalculation_(
  const double&, const int& , const int& , 
#if (NDIM>2)
  const int& , 
#endif
  const double*,const double*,
#if (NDIM>2)
  const double*, 
#endif
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*, 
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
#if (NDIM>1)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      ); 
 
#if (NDIM>2)
  void fluxcorrec2d_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const double*, const double*, const double *,const int&   ,
  const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 

  void fluxcorrec3d_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const double*, const double*, const double *,
  const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 
#endif

#if (NDIM==2)
  void fluxcorrec_(
  const double&, const int& , const int& , const int& ,const int& ,
  const double*, const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*      );
#endif

  void   consdiff_(
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*, 
  const double*, 
#if (NDIM>1)
  const double*,
#endif
#if (NDIM>2)
  const double*,
#endif
  double*      );

  void getbdry_( const int& ,
  const int& , const int& , const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , const int& , const int& , 
#endif
  const int& ,
#if (NDIM>1)
  const int& , 
#endif 
#if (NDIM>2)
  const int& , 
#endif 
  const int& ,
  const double*, const double&,
  double*      , 
  const double*, const double*, const int&);

#if (NDIM>2)
  void onethirdstate_(
  const double&, const double*, const double*, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, 
  const double*, const double*, const double*, 
  double*      ); 

  void fluxthird_(
  const double&, const double*, const double*, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, 
  const double*, 
  double*      , double*      , double*      );

  void fluxcorrecjt_(
  const double&, const double*, const double*, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, 
  const double*, const double*, const double*,
  double*      , double*      , double*      ,
  double*      , double*      , double*      );
#endif

   void detectgrad_(
#if (NDIM == 2)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
#if (NDIM == 3)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , 
      const int&, const int&,
      const double*,
      int* );

   void stufprobc_(
      const int& , const int& , const int&,
      const int& , const int& , const int& );
 
}
