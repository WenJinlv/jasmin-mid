define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      
          subroutine cjstabledt(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  gamma,cjsv,density,velocity,pressure,df,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &        ngc0,ngc1,ngc2
c
      REAL  
     &  gamma,cjsv,
     &  density(CELL3dVECG(ifirst,ilast,ngc)),
     &  velocity(CELL3dVECG(ifirst,ilast,ngc),0:NDIM-1),
     &  pressure(CELL3dVECG(ifirst,ilast,ngc)),
     &  df(CELL3dVECG(ifirst,ilast,ngc))
c    
      integer ic0,ic1,ic2,k

      REAL maxspeed(0:NDIM-1),lambda

      do k=0,NDIM-1
         maxspeed(k)=zero
      enddo

      stabdt=1.0d0

      do  ic2=ifirst2,ilast2
        do  ic1=ifirst1,ilast1
          do  ic0=ifirst0,ilast0
            if(df(ic0,ic1,ic2).ge.0.9998) then
               lambda = sqrt(max(zero,gamma*pressure(ic0,ic1,ic2)
     &                                      /density(ic0,ic1,ic2)))
               do k=0,NDIM-1
                  maxspeed(k)=max(maxspeed(k),
     &                abs(velocity(ic0,ic1,ic2,k))+lambda)
                  stabdt = min(stabdt,dx(k)/maxspeed(k))
               enddo
            else if(df(ic0,ic1,ic2).ge.1.d-4) then
               lambda = sqrt(max(zero,gamma*pressure(ic0,ic1,ic2)
     &                                      /density(ic0,ic1,ic2)))
               do k=0,NDIM-1
                  maxspeed(k)=max(maxspeed(k),
     &                abs(velocity(ic0,ic1,ic2,k))+lambda+5.0*cjsv)
                  stabdt = min(stabdt,dx(k)/maxspeed(k))
               enddo
            endif
          enddo
        enddo
      enddo

      return
      end 

