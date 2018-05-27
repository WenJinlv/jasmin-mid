define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine stabledt(dx,xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  uval,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL 
     &stabdt,dx(0:NDIM-1),xlo(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ngc0,ngc1
c
      REAL  
     &  uval(CELL2dVECG(ifirst,ilast,ngc))
      REAL 
     &  maxspeed(0:NDIM-1),vel(0:NDIM-1)
c
      integer ic0,ic1
c    
      maxspeed(0)=zero
      maxspeed(1)=zero
c   
      do ic1=ifirst1,ilast1
      vel(0)=-(xlo(1)+dx(1)*(ic1-ifirst1+half))
      do ic0=ifirst0,ilast0
         vel(1)=(xlo(0)+dx(0)*(ic0-ifirst0+half))
         maxspeed(0) = max(maxspeed(0), abs(vel(0)))
         maxspeed(1) = max(maxspeed(1), abs(vel(1)))
      enddo
      enddo
      stabdt=1.d+10
      do ic0 = 0,NDIM-1
         if(maxspeed(ic0).gt.1.d-6) then
            stabdt=min(stabdt,dx(ic0)/maxspeed(ic0))
         endif
      enddo
      return
      end 
