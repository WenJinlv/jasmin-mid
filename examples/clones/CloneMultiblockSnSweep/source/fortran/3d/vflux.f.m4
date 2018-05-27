define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim3d.i)dnl
c***********************************************************************
      subroutine vflux(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     mg0, ms0,
     &     wf, fi,
     &     sw,
     &     mg_begin, mg_end, mg_end2, loc2glob )
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
c
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2
      integer loc2glob(mg_begin:mg_end2), ig

      REAL
     &     wf    (NODE3d(ifirst,ilast,0), 0:mg0-1),
     &     fi    (mg_begin:mg_end2,NODE3d(ifirst,ilast,0),0:ms0-1)
      REAL 
     &     sw(0:ms0-1)
c    
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
c
      integer n, i,j,k, mg, ms
c
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg = 0, mg0-1
        wf(i,j,k,mg) = 0.0d0
      enddo  
      enddo  
      enddo  
      enddo  
c
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg = mg_begin, mg_end
         ig = loc2glob( mg )
         do ms = 0, ms0-1
            wf(i,j,k,ig)=wf(i,j,k,ig) + sw(ms)*fi(mg,i,j,k,ms)
         enddo
      enddo
      enddo
      enddo
      enddo
c
      return
      end   
c
c***********************************************************************
c
      subroutine vfluxerror(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     mg0,
     &     wf_0, wf_1, error,
     &     mg_begin, mg_end, mg_end2)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2
      REAL
     &     wf_0(NODE3d(ifirst,ilast,0), 0:mg0-1),
     &     wf_1(NODE3d(ifirst,ilast,0), 0:mg0-1)
      REAL error
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer n, i,j,k, mg, ms
      REAL a,b,ab,aberr

      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg= 0, mg0-1
         a = wf_0(i,j,k,mg)
         b = wf_1(i,j,k,mg)
         ab = max(abs(a),abs(b))
         if ( ab .ne. 0d0 ) aberr=(a-b)/ab
         if ( error .lt. aberr ) error = aberr
      enddo
      enddo
      enddo
      enddo
c
      return
      end

c***********************************************************************

