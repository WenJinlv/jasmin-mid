//
// 文件名:     EulerFort.h
// 软件包:     JASMIN application
// 版权  :     (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:     $Revision: 170 $
// 修改  :     $Date: 2007-06-27 08:44:22  $
// 描述  :     求解Euler方程组的F77子程序.
//

extern "C" {

  void eulerinit_(
  const int& , const double*, const double*, const double*,
  const int& , const int& ,
  const int& , const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int&,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,
  double*      , double*      , double*      ,
  const int&,
  const double*,
  const double*, const double*, const double*);
  
  void cjeulerinit_(
  const int& , const double*, const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int&,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  double*, double*, double*, double *,
  const double*, const double*, 
  const double&, const double&, const double&, const double&);


  void eulerinitsphere_(
  const int& , const double*, const double*, const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int&,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,
  double*      , double*      , double*      ,
  const double&, const double*, const double&,
  const double&, const double*, const double&,
  const double*, const double&);
  
  void   cjstabledt_(
  const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& ,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,const double&, 
  const double*, const double*, const double*, const double*,
  double&);
 


  void   cjconsdiff_(
  const int&, const double *, const double *,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*,
  const double*,
#if (NDIM>2)
  const double*,
#endif
   double*, double*, double*, double*,
  const double&, const double&, const double&, const double&, 
  const double&, const double&, const double&, 
  const double*, const double*, const double&); 

  void   consdiff_(
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*,
  const double*,
  const double*,
#if (NDIM>2)
  const double*,
#endif
  const double&, 
  double*      , double*      , double*      );



#if (NDIM == 2)
   void conservlinint2d_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*,
      double*, double*, double*, const int&,
      double*, double*, double*,
      double*, double*,
      double*, double*, double*, double*);

   void conservavg2d_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*);
#endif
#if (NDIM == 3)
   void conservlinint3d_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*,
      double*, double*, double*, const int&,
      double*, double*, double*,
      double*, double*, double*,
      double*, double*, double*, double*, double*, double*);

   void conservavg3d_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*);
#endif

   void cjdetectgrad_(
      const int& , const int& , 
      const int& , const int& , 
#if (NDIM>2)
      const int& , const int& , 
#endif
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#if (NDIM>2)
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , 
      const int&, const int&,
      const double*,
      const double*,
      const int* );
      
   void detectgrad_(
      const int& , const int& , 
      const int& , const int& , 
#if (NDIM>2)
      const int& , const int& , 
#endif
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#if (NDIM>2)
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void detectshock_(
      const int& , const int& , 
      const int& , const int& , 
#if (NDIM>2)
      const int& , const int& , 
#endif
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#if (NDIM>2)
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void stufprobc_(
      const int& , const int& , const int& , 
      const int& , const int& , const int& , 
      const int& , const int& , const int& ,
      const int& , const int& , const int& );

}
