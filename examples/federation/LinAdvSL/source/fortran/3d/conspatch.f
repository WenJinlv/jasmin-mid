


      subroutine conspatch(dx,xlo,
     &                     ifirst0,ilast0,
     &                     ifirst1,ilast1,
     &                     ifirst2,ilast2,
     &                     gwd0,gwd1,gwd2,
     &                     uold,flux0,flux1,flux2,unew)
      implicit none
      real*8 dx(0:2), xlo(0:2)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2

      real*8  uold(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2),
     &        flux0(ifirst0:ilast0+1,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2),
     &        flux1(ifirst0:ilast0,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2),
     &        flux2(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2+1),
     &        unew(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2)

      integer i,j,k

      do k=ifirst2,ilast2
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         unew(i,j,k)=uold(i,j,k)+(flux0(i,j,k)-flux0(i+1,j,k))/dx(0)
     &                          +(flux1(i,j,k)-flux1(i,j+1,k))/dx(1)
     &                          +(flux2(i,j,k)-flux2(i,j,k+1))/dx(2)
      enddo
      enddo
      enddo

      return
      end
               
