
c***********************************************************************
      subroutine vflux(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     mg0, ms0,
     &     wf, fi,
     &     sw)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      double precision
     &     wf    (ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:mg0-1),
     &     fi    (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1)
      double precision 
     &     sw(0:ms0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer n, i,j,k, mg, ms

      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg= 0, mg0-1
         wf(i,j,k,mg)=0d0
         do ms =0, ms0-1
            wf(i,j,k,mg)=wf(i,j,k,mg)+sw(ms)*fi(mg,i,j,k,ms)
         enddo
      enddo
      enddo
      enddo
      enddo

      return
      end   

c***********************************************************************

      subroutine vfluxerror(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     mg0,
     &     wf_0, wf_1, error)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      double precision
     &     wf_0(ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:mg0-1),
     &     wf_1(ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:mg0-1)
      double precision error
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer n, i,j,k, mg, ms
      double precision a,b,ab,aberr

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

      return
      end

c***********************************************************************

