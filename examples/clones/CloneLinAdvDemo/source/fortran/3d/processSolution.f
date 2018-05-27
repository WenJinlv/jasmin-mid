      subroutine prepsolution(ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        clone_number,
     &                        uval)
      implicit none
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2,clone_number

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  ratio

      ratio = clone_number+1
      if(clone_number.le.0) ratio=1.0d0
      do k=ifirst2-gwd2,ilast2+gwd2
         do j=ifirst1-gwd1,ilast1+gwd1
            do i=ifirst0-gwd0,ilast0+gwd0
               uval(i,j,k) = uval(i,j,k)*ratio
            enddo
         enddo
      enddo

      return
      end

      subroutine initpostsolution(ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        number_clones,
     &                        uval,uval_record)
      implicit none
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2,number_clones

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)
      real*8  uval_record(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  ratio_0,ratio_1

      ratio_0=1.0d0
      ratio_1=1.0d0
      if(number_clones>0) then
         ratio_0 = number_clones*(number_clones+1)/2.0d0
      endif

      do k=ifirst2-gwd2,ilast2+gwd2
         do j=ifirst1-gwd1,ilast1+gwd1
            do i=ifirst0-gwd0,ilast0+gwd0
               uval(i,j,k) = uval(i,j,k)/ratio_0
               uval_record(i,j,k) = uval_record(i,j,k)/ratio_1
            enddo
         enddo
      enddo

      return
      end

      subroutine postsolution(ifirst0,ilast0,
     &                        ifirst1,ilast1,
     &                        ifirst2,ilast2,
     &                        gwd0,gwd1,gwd2,
     &                        number_clones,
     &                        uval,uval_record)
      implicit none
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1,gwd2,number_clones

      real*8  uval(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)
      real*8  uval_record(ifirst0-gwd0:ilast0+gwd0,
     &             ifirst1-gwd1:ilast1+gwd1,
     &             ifirst2-gwd2:ilast2+gwd2)

      integer i,j,k
      real*8  ratio_0,ratio_1

      ratio_0=1.0d0
      ratio_1=1.0d0
      if(number_clones>0) then
         ratio_0 = number_clones*(number_clones+1)/2.0d0
         ratio_1 = number_clones*(number_clones+1)
     &                          *(2*number_clones+1)/6.0d0
      endif

      do k=ifirst2-gwd2,ilast2+gwd2
         do j=ifirst1-gwd1,ilast1+gwd1
            do i=ifirst0-gwd0,ilast0+gwd0
               uval(i,j,k) = uval(i,j,k)/ratio_0
               uval_record(i,j,k) = uval_record(i,j,k)/ratio_1
            enddo
         enddo
      enddo

      return
      end
