      subroutine initcoords0(
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        gwd0,gwd1,
     &                        coords)
      implicit none
      real*8 pi
      parameter(pi=3.1415926)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      real*8  coords(ifirst0-gwd0:ilast0+1+gwd0,
     &             ifirst1-gwd1:ilast1+1+gwd1,
     &             0:1)

      integer i,j

      do j=ifirst1,ilast1+1
      do i=ifirst0,ilast0+1
         coords(i,j,0)=i*1.6
         coords(i,j,1)=j*1.6
      enddo
      enddo

      return
      end
