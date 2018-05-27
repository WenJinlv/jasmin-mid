      subroutine initsetvalue(dx,xlo,
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        uval)
      implicit none
      real*8 pi
      parameter(pi=3.1415926d0)
      real*8 dx(0:2), xlo(0:2)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  xhalf, yhalf, zhalf

      do k=ifirst2,ilast2
         zhalf = xlo(2)+(k-ifirst2+0.5d0)*dx(2)
         do j=ifirst1,ilast1
            yhalf = xlo(1)+(j-ifirst1+0.5d0)*dx(1)
            do i=ifirst0,ilast0
               xhalf = xlo(0)+(i-ifirst0+0.5d0)*dx(0)
               if(xhalf.le.(2.0-dsin(pi*yhalf))) then 
                  uval(i,j,k)=2.0d0
               else
                  uval(i,j,k)=0.0d0
               endif
            enddo
         enddo
      enddo

      return
      end
