
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine tagpatch(dx,xlo,
     &                    ifirst0,ilast0,
     &                    ifirst1,ilast1,
     &                    gwd0,gwd1,
     &                    uval,tag)
      implicit none
      real*8 dx(0:1), xlo(0:1)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      real*8  uval(CELL2dVECG(ifirst,ilast,gwd))
      integer tag(CELL2d(ifirst,ilast,0))

      integer i,j
      real*8  gradl, gradr, grad, threshold

      threshold = 1.d0
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         gradl = abs(uval(i,j)-uval(i-1,j))/dx(0)
         gradr = abs(uval(i+1,j)-uval(i,j))/dx(0)
         grad = dmax1(gradl,gradr)
         tag(i,j)=0
         if(grad.gt.threshold) tag(i,j)=1
      enddo
      enddo

      return
      end
