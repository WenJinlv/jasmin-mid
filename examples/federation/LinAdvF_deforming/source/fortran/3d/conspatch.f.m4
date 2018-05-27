
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine conspatch(
     &                     ifirst0,ilast0,
     &                     ifirst1,ilast1,
     &                     ifirst2,ilast2,
     &                     gwd0,gwd1,gwd2,
     &                     uold,flux0,flux1,flux2,unew)
      implicit none
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2

      real*8  uold(CELL3d(ifirst,ilast,0)),
     &        flux0(SIDE3d0(ifirst,ilast,0)),
     &        flux1(SIDE3d1(ifirst,ilast,0)),
     &        flux2(SIDE3d2(ifirst,ilast,0)),
     &        unew(CELL3d(ifirst,ilast,0))

      integer i,j,k

      do k=ifirst2,ilast2
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         unew(i,j,k)=uold(i,j,k)+(flux0(i,j,k)-flux0(i+1,j,k))
     &                          +(flux1(i,j,k)-flux1(i,j+1,k))
     &                          +(flux2(i,j,k)-flux2(i,j,k+1))
      enddo
      enddo
      enddo

      return
      end
               
