//
// 文件名: RotAdvFort.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.2 $
// 修改  : $Date: 2007/01/06 09:07:25 $
// 描述  : F77外部子程序接口.
//


#if (NDIM==2)
#define BOX_CLAIM                                   \
        const int &, const int &,           \
        const int &, const int &
#define BOX_CLAIM1                                   \
        const int &ifirst01, const int &ilast01,           \
        const int &ifirst11, const int &ilast11
#define ABC_CLAIM                                \
        double *, double *
#elif (NDIM==3)
#define BOX_CLAIM                                   \
        const int &, const int &,           \
        const int &, const int &,           \
        const int &, const int &
#define BOX_CLAIM1                                   \
        const int &ifirst01, const int &ilast01,           \
        const int &ifirst11, const int &ilast11,           \
        const int &ifirst21, const int &ilast21
#define ABC_CLAIM                                \
        double *, double *, double *
#endif

extern "C" {
 
   void snquad_(double *angles,
                    const double *sw,
                    const double *u0,
                    const int &ms0, const int &n0);

   void init_(BOX_CLAIM,
               const int &mg0, const int &mg1, const int &ms0,
               double *sgm_t, double *sgm_s, double *wf, double *fi,
               const double *dx, const double *xleng,
               const double *sgm_t0, const double *sgm_s0,
               const double *vn,
               const double *angles, const double *sw,
               ABC_CLAIM);

   void scatter_(double *qs,
               const double *sgm_s,
               const double *wf,
               BOX_CLAIM,
               const int &mg0,   const int &mg1);
               

   void sourcew_(BOX_CLAIM,
               const int &mg0, const int &ms0,
               double *qw,
               const double *dx, const double *xleng,
               const double *angles,
               const double *vn,
               const double *ysgm_s, const double *ysgm_t, 
               const double &yd,
               ABC_CLAIM);
               
   void sourceqnw_(BOX_CLAIM,
               const int &mg0, const int &ms0,
               const double &dt,
               double *qnw,
               const double *qw, const double *fi,
               const double *vn);
               
   void solve_(BOX_CLAIM,
               const int &mg0, const int &ms0,
               const int &ms,  const double *angle,
               const int &ncells, const int cells[][NDIM],
               double *fi,
               const double *qnw, const double *qs, 
               const double *sgm_t,
               const double &dt, const double *dx,
               const double *vn ); 
    
  void setdependency_(BOX_CLAIM,
               const int &ms0,
               int *depend0, int *depend1,
#if (NDIM==3)
               int *depend2,
#endif
               const double *angles);

   void vflux_(BOX_CLAIM,
               const int &mg0, const int &ms0,
               double *wf,
               const double *fi,
               const double *sw);

   void vfluxerror_(BOX_CLAIM,
               const int &mg0, 
               const double *wf_0, const double *wf_1, double *error);
 
   void computeenergy_( BOX_CLAIM, BOX_CLAIM1,
                const int &mg0, 
                const double *wf,
                double *buffer );

   void diffsolution_(BOX_CLAIM,
               const int &mg0, const int &ms0,
               const double *fi,
               const double *dx, const double *xleng,
               const double &yd, const double *vn, const double *angles,
               ABC_CLAIM);
  
}
