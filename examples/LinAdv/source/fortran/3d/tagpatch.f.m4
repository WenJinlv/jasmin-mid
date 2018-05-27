
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine tagpatch(dx,xlo,
     &                    ifirst0,ilast0,
     &                    ifirst1,ilast1,
     &                    ifirst2,ilast2,
     &                    gwd0,gwd1,gwd2,
     &                    uval,tag)
      implicit none
      real*8 dx(0:2), xlo(0:2)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1, gwd2

      real*8  uval(CELL3dVECG(ifirst,ilast,gwd))
      integer tag(CELL3d(ifirst,ilast,0))

      integer i,j,k
      real*8  gradl, gradr, grad, threshold

      threshold = 1.d0
      do k=ifirst2,ilast2
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         gradl = abs(uval(i,j,k)-uval(i-1,j,k))/dx(0)
         gradr = abs(uval(i+1,j,k)-uval(i,j,k))/dx(0)
         grad = dmax1(gradl,gradr)
         tag(i,j,k)=0
         if(grad.gt.threshold) tag(i,j,k)=1
      enddo
      enddo
      enddo

      return
      end
