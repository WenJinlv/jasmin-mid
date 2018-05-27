
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine advancepatch(velocity,dt,
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        uval,flux0,flux1,flux2)
      implicit none
      real*8 velocity,dt
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2

      real*8  uval(CELL3dVECG(ifirst,ilast,gwd)),
     &        flux0(SIDE3d0(ifirst,ilast,0)),
     &        flux1(SIDE3d1(ifirst,ilast,0)),
     &        flux2(SIDE3d2(ifirst,ilast,0))

      integer i,j,k

      do k=ifirst2,ilast2
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0+1
         flux0(i,j,k)=dt*velocity*uval(i-1,j,k)
      enddo
      enddo
      enddo

      do k=ifirst2,ilast2
      do j=ifirst1,ilast1+1
      do i=ifirst0,ilast0
         flux1(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=ifirst2,ilast2+1
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         flux2(i,j,k)=0.0
      enddo
      enddo
      enddo

      return
      end
               
