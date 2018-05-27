      subroutine initsetvalue(dx,xlo,cnid,
     &                        ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        uval,uval_record)
      implicit none
      real*8 pi
      parameter(pi=3.1415926)
      real*8 dx(0:2), xlo(0:2)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2,cnid

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)
      real*8  uval_record(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  xhalf, yhalf, zhalf,ratio

      ratio=cnid*1.0d0
      if(cnid.le.0) ratio=0.0d0
      do k=ifirst2-gwd2,ilast2+gwd2
         zhalf = xlo(2)+(k-ifirst2+0.5d0)*dx(2)
         do j=ifirst1-gwd1,ilast1+gwd1
            yhalf = xlo(1)+(j-ifirst1+0.5d0)*dx(1)
            do i=ifirst0-gwd0,ilast0+gwd0
               xhalf = xlo(0)+(i-ifirst0+0.5d0)*dx(0)
               if(xhalf.le.(2.0d0-dsin(pi*yhalf))) then 
                  uval(i,j,k)=2.0d0*(ratio+1)
               else
                  uval(i,j,k)=0.0d0
               endif
            enddo
         enddo
      enddo

      do k=ifirst2-gwd2,ilast2+gwd2
      do j=ifirst1-gwd1,ilast1+gwd1
      do i=ifirst0-gwd0,ilast0+gwd0
         uval_record(i,j,k)=uval(i,j,k)
      enddo
      enddo
      enddo

      return
      end
               
