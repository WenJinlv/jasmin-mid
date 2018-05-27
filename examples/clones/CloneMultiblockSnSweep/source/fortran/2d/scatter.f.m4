define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim2d.i)dnl
c***********************************************************************
      subroutine scatter(
     &  qs,sgm_s,wf,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0,mg1, mg_begin, mg_end, mg_end2,
     &  loc2glob )
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0,mg1
      REAL
     &     qs   (mg_begin:mg_end2,NODE2d(ifirst,ilast,0)),
     &     sgm_s(0:mg1-1,NODE2d(ifirst,ilast,0)),
     &     wf   ( NODE2d(ifirst,ilast,0), 0:mg0-1 )
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j, loc2glob(mg_begin:mg_end2)
      integer mg, ig, lg1
      integer mg_begin, mg_end, mg_end2
      integer k, kg_glob, ig_begin, mg_glob
c
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         lg1 = 0
         ig_begin = 0
         do mg = mg_begin, mg_end
            mg_glob = loc2glob( mg ) 
            do kg_glob = ig_begin, mg_glob - 1
               do k=0, kg_glob
                  lg1 = lg1+1
               enddo
            enddo

            qs(mg,i,j)=0d0
            do ig=0, mg_glob
               qs(mg,i,j)=qs(mg,i,j)+sgm_s(lg1,i,j)*wf(i,j,ig)
               lg1=lg1+1
            enddo
            ig_begin = mg_glob+1

         enddo
      enddo
      enddo
c
      return
      end   

c***********************************************************************
