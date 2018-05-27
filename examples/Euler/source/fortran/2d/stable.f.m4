define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      
      subroutine cjstabledt(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  gamma,cjsv,density,velocity,pressure,df,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ngc0,ngc1
c
      REAL  
     &  gamma,cjsv,
     &  density(CELL2dVECG(ifirst,ilast,ngc)),
     &  velocity(CELL2dVECG(ifirst,ilast,ngc),0:NDIM-1),
     &  pressure(CELL2dVECG(ifirst,ilast,ngc)),
     &  df(CELL2dVECG(ifirst,ilast,ngc))
c    
      integer ic0,ic1

      REAL maxspeed(0:NDIM-1),lambda

      maxspeed(0)=zero
      maxspeed(1)=zero

      stabdt=1.0d0

      do  ic1=ifirst1,ilast1
         do  ic0=ifirst0,ilast0
            if(df(ic0,ic1).ge.0.9998) then
               lambda = 
     &          sqrt(max(zero,gamma*pressure(ic0,ic1)/density(ic0,ic1)))
               maxspeed(0) = max(maxspeed(0),
     &              abs(velocity(ic0,ic1,0))+lambda)
               maxspeed(1) = max(maxspeed(1),
     &              abs(velocity(ic0,ic1,1))+lambda)
               stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
            else if(df(ic0,ic1).ge.1.d-4) then
               lambda =
     &          sqrt(max(zero,gamma*pressure(ic0,ic1)/density(ic0,ic1)))
               maxspeed(0) = max(maxspeed(0),
     &              abs(velocity(ic0,ic1,0))+lambda+5.0*cjsv)
               maxspeed(1) = max(maxspeed(1),
     &              abs(velocity(ic0,ic1,1))+lambda+5.0*cjsv)
               stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
            endif
         enddo
      enddo

      return
      end       
