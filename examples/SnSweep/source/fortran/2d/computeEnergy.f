









c***********************************************************************
      subroutine computeenergy(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     ibeg0,iend0,ibeg1,iend1,
     &     mg0, 
     &     wf, 
     &     energy)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer mg0
      double precision
     &     wf    (ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:mg0-1),
     &     energy(ibeg0:iend0+1,
     &          ibeg1:iend1+1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j, mg

      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         energy(i,j) = 0d0
         do mg= 0, mg0-1
            energy(i,j) = energy(i,j) + wf(i,j,mg)
         enddo
      enddo
      enddo

      return
      end   

c***********************************************************************

