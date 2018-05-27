define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim2d.i)dnl
c***********************************************************************
      subroutine scatter(
     &  qs,sgm_s,wf,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0,mg1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0,mg1
      REAL
     &     qs   (0:mg0-1,NODE2d(ifirst,ilast,0)),
     &     sgm_s(0:mg1-1,NODE2d(ifirst,ilast,0)),
     &     wf   (NODE2d(ifirst,ilast,0),0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j
      integer mg, ig, lg1
     
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         lg1=0
         do mg=0,mg0-1
            qs(mg,i,j)=0d0
            do ig=0,mg
               qs(mg,i,j)=qs(mg,i,j)+sgm_s(lg1,i,j)*wf(i,j,ig)
               lg1=lg1+1
            enddo
         enddo
      enddo
      enddo

      return
      end   

c***********************************************************************
