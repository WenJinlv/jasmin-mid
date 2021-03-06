      subroutine pulseclonecomputing(clone_number,
     &                     src_clone_number,
     &                     ifirst0,ilast0,
     &                     ifirst1,ilast1,
     &                     ifirst2,ilast2,
     &                     gwd0,gwd1,gwd2,
     &                     uval,uval_record)
      implicit none
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2
      integer clone_number, src_clone_number

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)
      real*8  uval_record(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  ratio

      if(clone_number.le.0) then

         if(src_clone_number.le.0) then
            do k=ifirst2-gwd2,ilast2+gwd2
            do j=ifirst1-gwd1,ilast1+gwd1
            do i=ifirst0-gwd0,ilast0+gwd0
               uval_record(i,j,k)=0.0d0
            enddo
            enddo
            enddo
         endif

         ratio = src_clone_number+1
         if(src_clone_number.le.0) ratio = 1.0d0
         do k=ifirst2-gwd2,ilast2+gwd2
         do j=ifirst1-gwd1,ilast1+gwd1
         do i=ifirst0-gwd0,ilast0+gwd0
            uval_record(i,j,k)=uval_record(i,j,k)+uval(i,j,k)*ratio
         enddo
         enddo
         enddo

      endif

      return
      end
               
