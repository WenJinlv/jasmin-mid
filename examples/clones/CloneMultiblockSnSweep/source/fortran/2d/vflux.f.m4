define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim2d.i)dnl
c***********************************************************************
      subroutine vflux(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     mg0, ms0,
     &     wf, fi,
     &     sw,
     &     mg_begin, mg_end, mg_end2, loc2glob )
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
c
      integer ifirst0, ilast0, ifirst1, ilast1
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2, ig
      integer loc2glob(mg_begin:mg_end2) 
c
      REAL
     &     wf    (NODE2d(ifirst,ilast,0),0:mg0-1),
     &     fi    (mg_begin:mg_end2,NODE2d(ifirst,ilast,0), 0:ms0-1)
      REAL 
     &     sw(0:ms0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer n, i,j, mg, ms
c
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg= 0, mg0-1
         wf(i,j,mg)=0d0
      enddo
      enddo
      enddo
c
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg= mg_begin, mg_end
        ig = loc2glob( mg )
        do ms = 0, ms0-1
          wf(i,j,ig)=wf(i,j,ig)+sw(ms)*fi(mg,i,j,ms)
        enddo
      enddo
      enddo
      enddo
c      
      return
      end   
c***********************************************************************

      subroutine vfluxerror(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     mg0,
     &     wf_0, wf_1, error,
     &     mg_begin, mg_end, mg_end2 )
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2 
      REAL
     &     wf_0( NODE2d(ifirst,ilast,0), 0:mg0-1 ),
     &     wf_1( NODE2d(ifirst,ilast,0), 0:mg0-1 )
      REAL error
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer n, i,j, mg, ms
      REAL a,b,ab, aberr

      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg= 0,mg0-1
         a = wf_0(i,j,mg)
         b = wf_1(i,j,mg)
         ab = max(abs(a),abs(b))
         if ( ab .ne. 0d0 ) aberr=(a-b)/ab
         if ( error .lt. aberr ) error = aberr
      enddo
      enddo
      enddo
c
      return
      end

c***********************************************************************

