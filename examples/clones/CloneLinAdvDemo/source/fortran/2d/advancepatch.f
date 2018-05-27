      subroutine advancepatch(velocity,dt,dx,xlo,
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        gwd0,gwd1,
     &                        uval,flux0,flux1)
      implicit none
      real*8 dx(0:1), xlo(0:1),velocity,dt
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1),
     &        flux0(ifirst0:ilast0+1, ifirst1:ilast1),
     &        flux1(ifirst0:ilast0  , ifirst1:ilast1+1)

      integer i,j

      do j=ifirst1,ilast1
      do i=ifirst0,ilast0+1
         flux0(i,j)=dt*velocity*uval(i-1,j)
      enddo
      enddo

      do j=ifirst1,ilast1+1
      do i=ifirst0,ilast0
         flux1(i,j)=0.0
      enddo
      enddo

      return
      end
               
