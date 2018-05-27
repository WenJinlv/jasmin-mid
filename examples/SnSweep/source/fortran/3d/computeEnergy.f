
c***********************************************************************
      subroutine computeenergy(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &     mg0, 
     &     wf, 
     &     energy)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer mg0
      double precision
     &     wf    (ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:mg0-1),
     &     energy(ibeg0:iend0+1,
     &          ibeg1:iend1+1,
     &          ibeg2:iend2+1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j,k, mg

      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         energy(i,j,k) = 0d0
         do mg= 0, mg0-1
            energy(i,j,k) = energy(i,j,k) + wf(i,j,k,mg)
         enddo
      enddo
      enddo
      enddo

      return
      end   

c***********************************************************************

