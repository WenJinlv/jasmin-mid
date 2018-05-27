//
// 文件名:  MolecularDynamicsFort.h
// 软件包:  JASMIN application
// 版权  :  (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号:  $Revision: 170 $
// 修改  :  $Date: 2007-06-27 08:44:22 +0800 (涓, 27  6 2007) $
// 描述  :  F77数值计算子程序接口.
//

extern "C" {

    void setglobalparam_(const int*, const int*,
                         const double*, const double*, const double*,
                         const int& , const int& ,
                         const int& , const int& ,
#if (NDIM>2)
                         const int& , const int& ,
#endif
                         const double*, const double*,
                         const double*, const double*,
                         const double*, const double*
                        );

    void parinit_(
        const double*, const double*, const double*,
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
        const int*, const int* ,
        const double*, const int*, 
        const int*, const int*, 
        const int*,

        const int&,
        const int& ,
#if (NDIM>2)
        const int& ,
#endif
        const double*, const int*
    );


    void initialload_(
        const double*, const double*, const double*,
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
        const double*
    );

    void computeload_(
        const double*, const double*, 
        const double*, const double*,
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
        const int&,
        const int& ,
#if (NDIM>2)
        const int& ,
#endif
        const int*, const double*
    );

    void computedensity_(
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
        const int*, 
        const int&,
        const int& ,
#if (NDIM>2)
        const int& ,
#endif
       const double*
    );

    void updateparticles_(
        const double* , const double*, 
        const double*, const double*,
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
        const int*,  const int*, 
        const double*, const int*,
        const int*, const int*, const int*, 
        double *
    );

}
