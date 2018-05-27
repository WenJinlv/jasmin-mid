define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine stabledt(omega,dx,xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  uval,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,omega(0:NDIM-1),dx(0:NDIM-1),xlo(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     & ngc0,ngc1,ngc2
c
      REAL
     &  uval(CELL3dVECG(ifirst,ilast,ngc))
c
      integer ic0,ic1,ic2
      REAL maxspeed(0:NDIM-1),vel(0:NDIM-1)
c
      maxspeed(0)=zero
      maxspeed(1)=zero
      maxspeed(2)=zero

      do ic2=ifirst2,ilast2
      do ic1=ifirst1,ilast1
      do ic0=ifirst0,ilast0
         vel(0)= (xlo(1)+dx(1)*(ic1-ifirst1+half))*omega(2)
     &          -(xlo(2)+dx(2)*(ic2-ifirst2+half))*omega(1)
         vel(1)= (xlo(2)+dx(2)*(ic2-ifirst2+half))*omega(0)
     &          -(xlo(0)+dx(0)*(ic0-ifirst0+half))*omega(2)
         vel(2)= (xlo(0)+dx(0)*(ic0-ifirst0+half))*omega(1)
     &          -(xlo(1)+dx(1)*(ic1-ifirst1+half))*omega(0)
         maxspeed(0) = max(maxspeed(0), abs(vel(0)))
         maxspeed(1) = max(maxspeed(1), abs(vel(1)))
         maxspeed(2) = max(maxspeed(2), abs(vel(2)))
      enddo
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
