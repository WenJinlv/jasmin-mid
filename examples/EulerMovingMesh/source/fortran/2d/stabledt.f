










c***********************************************************************
c
c    求稳定步长.
c
c***********************************************************************
      subroutine stabledt(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,gamma,
     &  xo, yo, uval,
     &  dt)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      double precision
     &    dt,gamma
c
c variables in 2d indexed
c
      double precision
     &     xo(ifirst0-gcw0:ilast0+1+gcw0,
     &          ifirst1-gcw1:ilast1+1+gcw1),
     &     yo(ifirst0-gcw0:ilast0+1+gcw0,
     &          ifirst1-gcw1:ilast1+1+gcw1),
     &     uval(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,2+2)
c
      integer i,j
      double precision
     &  sum, dx, dy, rho, pre, ccc, area, dtx, dty, alengx, alengy,
     &  aux, auy
c
c***********************************************************************     
c
      SUM=0.0
      DX=1.d20
      DY=1.d20    
      DO j=ifirst1,ilast1      
      DO i=ifirst0,ilast0
          rho=uval(i,j,1)
          aux=uval(i,j,2)/rho
          auy=uval(i,j,3)/rho
          PRE=(gamma-1.d0)*(uval(i,j,4)-0.5d0*rho*(aux*aux+auy*auy))
          CCC=SQRT(gamma*PRE/rho)
          SUM=MAX(SUM,ABS(aux)+CCC,ABS(auy)+CCC) 
          
          area=0.5*((xo(i,j)-xo(i+1,j+1))*(yo(i+1,j)-yo(i,j+1))
     &             -(yo(i,j)-yo(i+1,j+1))*(xo(i+1,j)-xo(i,j+1)))
          
          DTX=xo(i+1,j+1)-xo(i+1,j  )
          DTY=yo(i+1,j+1)-yo(i+1,j  )
          ALENGX=SQRT(DTX**2+DTY**2)
          DTX=xo(i  ,j+1)-xo(i+1,j+1)
          DTY=yo(i  ,j+1)-yo(i+1,j+1) 
          ALENGY=SQRT(DTX**2+DTY**2)                    
               
          DX=MIN(DX,area/ALENGX)
          DY=MIN(DY,area/ALENGY)
      ENDDO
      ENDDO
      
      DT=MIN(DX,DY)/SUM          
c
      return
      end
