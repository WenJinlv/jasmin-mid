









c***********************************************************************
      subroutine solve(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     mg0, ms0,
     &     ms, angles,
     &     ncells, cells,
     &     fi,
     &     qnw, qs, sgm_t,
     &     dt, dx, vn) 
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, ms0
      integer ms, ncells
      integer cells(0:2-1,ncells)
      double precision
     &     fi   (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1),
     &     qnw  (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1),
     &     qs   (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1),
     &     sgm_t(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1)
      double precision 
     &     dt, dx(0:2-1), 
     &     angles(0:2-1,0:ms0-1),
     &     vn(0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer n, i,j, mg, istep,jstep, i_up,j_up
      double precision    c1,c2,c0
     
      istep = 1
      jstep = 1
      if(angles(0,ms).lt.0d0) istep=0
      if(angles(1,ms).lt.0d0) jstep=0

      c1=dabs(angles(0,ms))/dx(0)
      c2=dabs(angles(1,ms))/dx(1)

      do n = 1, ncells
         i=cells(0,n)+istep
         j=cells(1,n)+jstep
         i_up = cells(0,n)+1-istep
         j_up = cells(1,n)+1-jstep
         do mg=0,mg0-1
            c0=1d0/(vn(mg)*dt)+c1+c2
            fi(mg,i,j,ms)=(qnw(mg,i,j,ms)
     1                      +c1*fi(mg,i_up,j,ms)
     2                      +c2*fi(mg,i,j_up,ms)
     4                      +qs(mg,i,j))/(c0+sgm_t(mg,i,j))
         enddo
      enddo

      return
      end   

c***********************************************************************

      subroutine setdependency(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ms0,
     &  depend0, depend1,
     &  angles)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ms0
      integer 
     &     depend0(ifirst0:ilast0+1,
     &          ifirst1:ilast1,0:ms0-1),
     &     depend1(ifirst0:ilast0,
     &          ifirst1:ilast1+1,0:ms0-1)
      double precision 
     &     angles(0:2-1,0:ms0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i0, i1, d0, d1, d
     
      do d = 0, ms0-1
         d0 = 0 
         if ( angles(0,d) .gt. 0d0 ) then
            d0 = 1
         else if ( angles(0,d) .lt. 0d0 ) then 
            d0 = -1
         endif

         d1 = 0 
         if ( angles(1,d) .gt. 0d0 ) then
            d1 = 1
         else if ( angles(1,d) .lt. 0d0 ) then 
            d1 = -1
         endif

         do i1 = ifirst1, ilast1
         do i0 = ifirst0, ilast0+1
            depend0(i0,i1,d) = d0
         enddo
         enddo

         do i1 = ifirst1, ilast1+1
         do i0 = ifirst0, ilast0
            depend1(i0,i1,d) = d1
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************
