


      subroutine advancepatch(velocity,dt,dx,xlo,
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        uval,flux0,flux1,flux2)
      implicit none
      real*8 dx(0:2), xlo(0:2),velocity,dt
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &          ifirst1-gwd1:ilast1+gwd1,
     &          ifirst2-gwd2:ilast2+gwd2),
     &        flux0(ifirst0:ilast0+1,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2),
     &        flux1(ifirst0:ilast0,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2),
     &        flux2(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2+1)

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
               
