define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim3d.i)dnl
c***********************************************************************
      subroutine scatter(
     &  qs,sgm_s,wf,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0,mg1, mg_begin, mg_end, mg_end2, loc2glob)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0,mg1
      REAL
     &     qs   (mg_begin:mg_end2,NODE3d(ifirst,ilast,0)),
     &     sgm_s(0:mg1-1, NODE3d(ifirst,ilast,0)),
     &     wf   (NODE3d(ifirst,ilast,0),0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j,k
      integer mg, ig, lg1
      integer mg_begin, mg_end, mg_end2
      integer loc2glob( mg_begin: mg_end2 )
      integer  n,ng_glob, ig_begin, mg_glob

     
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         lg1 = 0
         ig_begin = 0
         do mg = mg_begin, mg_end
            mg_glob = loc2glob( mg )
            do ng_glob = ig_begin, mg_glob - 1
               do n=0, ng_glob
                  lg1 = lg1+1
               enddo
            enddo
            qs(mg,i,j,k)=0d0
            do ig=0,mg_glob
              qs(mg,i,j,k)=qs(mg,i,j,k)+sgm_s(lg1,i,j,k)*wf(i,j,k,ig)
              lg1=lg1+1
            enddo
            ig_begin = mg_glob+1
         enddo
      enddo
      enddo
      enddo
c
      return
      end   
c***********************************************************************
