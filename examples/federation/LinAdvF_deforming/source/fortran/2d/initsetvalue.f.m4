      subroutine initsetvalue(
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        gwd0,gwd1,
     &                        uval)
      implicit none
      real*8 pi
      parameter(pi=3.1415926)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1)

      integer i,j

      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         uval(i,j)=i+8*j
      enddo
      enddo

      return
      end
