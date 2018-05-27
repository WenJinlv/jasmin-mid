
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine conspatch(dx,xlo,
     &                     ifirst0,ilast0,
     &                     ifirst1,ilast1,
     &                     gwd0,gwd1,
     &                     uold,flux0,flux1,unew)
      implicit none
      real*8 dx(0:1), xlo(0:1)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      real*8  uold(CELL2d(ifirst,ilast,0)),
     &        flux0(SIDE2d0(ifirst,ilast,0)),
     &        flux1(SIDE2d1(ifirst,ilast,0)),
     &        unew(CELL2d(ifirst,ilast,0))

      integer i,j

      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         unew(i,j)=uold(i,j)+(flux0(i,j)-flux0(i+1,j))/dx(0)
     &                       +(flux1(i,j)-flux1(i,j+1))/dx(1)
      enddo
      enddo

      return
      end
               
