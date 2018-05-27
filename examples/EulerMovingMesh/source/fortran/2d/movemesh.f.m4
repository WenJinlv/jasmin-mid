define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    移动网格.
c
c***********************************************************************
      subroutine movemesh(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  xo, yo, wo, xn, yn)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
c
c variables in 2d indexed
c
      REAL
     &     xo(NODE2dVECG(ifirst,ilast,gcw)),
     &     yo(NODE2dVECG(ifirst,ilast,gcw)),
     &     wo(CELL2dVECG(ifirst,ilast,gcw)),
     &     xn(NODE2d(ifirst,ilast,0)),
     &     yn(NODE2d(ifirst,ilast,0))
c
      integer J,K
      REAL
     &  CXP,CXM,CYP,CYM
c
c***********************************************************************     
C
      DO K=ifirst1,ilast1+1
      DO J=ifirst0,ilast0+1
         CXP=0.5*(WO(J,K)+WO(J,K-1))
         CXM=0.5*(WO(J-1,K)+WO(J-1,K-1))                       
         CYP=0.5*(WO(J,K)+WO(J-1,K))          
         CYM=0.5*(WO(J,K-1)+WO(J-1,K-1))
         XN(J,K)=( CXP*XO(J+1,K)+CXM*XO(J-1,K)
     &            +CYP*XO(J,K+1)+CYM*XO(J,K-1) )/(CXP+CXM+CYP+CYM)        
         YN(J,K)=( CXP*YO(J+1,K)+CXM*YO(J-1,K)
     &            +CYP*YO(J,K+1)+CYM*YO(J,K-1) )/(CXP+CXM+CYP+CYM)
      ENDDO
      ENDDO
c
      return
      end
