define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    计算控制函数.
c
c***********************************************************************
      subroutine computewcnew(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,alpha,beta,dxi,det,
     &  xo, yo, uo, wo)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL
     &    alpha,beta,dxi,det
c
c variables in 2d indexed
c
      REAL
     &     xo(NODE2dVECG(ifirst,ilast,gcw)),
     &     yo(NODE2dVECG(ifirst,ilast,gcw)),
     &     uo(CELL2dVECG(ifirst,ilast,gcw),NDIM+2),
     &     wo(CELL2d(ifirst,ilast,0))
      REAL
     &     FMonitor
c
      integer J,K
      REAL
     &  XXI,YXI,XET,YET,AREA,XIX,XIY,ETX,ETY,UU,UXI,UET,UX,UY
c
c***********************************************************************     
C
      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0
         XXI = 0.5*(XO(J+1,K)-XO(J-1,K))      
         YXI = 0.5*(YO(J+1,K)-YO(J-1,K)) 
         XET = 0.5*(XO(J,K+1)-XO(J,K-1))      
         YET = 0.5*(YO(J,K+1)-YO(J,K-1))              
         AREA= XXI*YET-XET*YXI
         XIX = YET/AREA
         XIY =-XET/AREA
         ETX =-YXI/AREA
         ETY = XXI/AREA
         UU=UO(J,K,1)
         UXI=(UO(J+1,K,1)-UO(J-1,K,1))/2.0/DXI
         UET=(UO(J,K+1,1)-UO(J,K-1,1))/2.0/DET   
         UX=UXI*XIX+UET*ETX
         UY=UXI*XIY+UET*ETY          
         WO(J,K)=FMonitor(UU,UXI,UET,ALPHA,BETA)   
      ENDDO
      ENDDO
C
      RETURN
      END
c
      FUNCTION FMonitor(X,DX,DY,ALPHA,BETA)
      IMPLICIT NONE
      REAL
     & X,DX,DY,ALPHA,BETA,FMonitor
C
C     DEFINE MONITOR FUNCTION
C
      FMonitor=SQRT(1.0+ALPHA*X*X+BETA*(DX*DX+DY*DY) )

      RETURN
      END
