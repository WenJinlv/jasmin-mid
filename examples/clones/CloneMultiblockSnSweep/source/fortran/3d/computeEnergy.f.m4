define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim3d.i)dnl
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
      REAL
     &     wf    (NODE3d(ifirst,ilast,0),0:mg0-1),
     &     energy(NODE3d(ibeg,iend,0))
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

