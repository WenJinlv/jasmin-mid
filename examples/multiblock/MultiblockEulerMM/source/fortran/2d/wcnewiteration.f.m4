define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c     迭代求解控制函数.
c
c***********************************************************************
      subroutine wcnewiteration(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,wo,wn)
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
     &     wo(CELL2dVECG(ifirst,ilast,gcw)),
     &     wn(CELL2d(ifirst,ilast,0))
c
      integer J,K
c
c***********************************************************************     
C
      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0
         WN(J,K)=WO(J,K)/4.0
     &        +(WO(J+1,K)+WO(J-1,K)+WO(J,K+1)+WO(J,K-1))/8.
     &        +(WO(J+1,K+1)+WO(J+1,K-1)+WO(J-1,K+1)+WO(J-1,K-1))/16.
      ENDDO
      ENDDO
c
      return
      end
