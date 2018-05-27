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

      real*8  uold (ifirst0-gwd0:ilast0+gwd0,
     &              ifirst1-gwd1:ilast1+gwd1),
     &        flux0(ifirst0:ilast0+1, ifirst1:ilast1),
     &        flux1(ifirst0:ilast0, ifirst1:ilast1+1),
     &        unew (ifirst0:ilast0, ifirst1:ilast1)

      integer i,j

      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         unew(i,j)=uold(i,j)+(flux0(i,j)-flux0(i+1,j))/dx(0)
     &                       +(flux1(i,j)-flux1(i,j+1))/dx(1)
      enddo
      enddo

      return
      end
               
