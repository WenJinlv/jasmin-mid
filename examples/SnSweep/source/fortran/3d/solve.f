
c***********************************************************************
      subroutine solve(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     mg0, ms0,
     &     ms, angles,
     &     ncells, cells,
     &     fi,
     &     qnw, qs, sgm_t,
     &     dt, dx, vn) 
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer ms, ncells
      integer cells(0:3-1,ncells)
      double precision
     &     fi   (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1),
     &     qnw  (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1),
     &     qs   (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1),
     &     sgm_t(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1)
      double precision 
     &     dt, dx(0:3-1), 
     &     angles(0:3-1,0:ms0-1),
     &     vn(0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer n, i,j,k, mg, istep,jstep,kstep, i_up,j_up,k_up
      double precision    c1,c2,c3,c0
     
      istep = 1
      jstep = 1
      kstep = 1
      if(angles(0,ms).lt.0d0) istep=0
      if(angles(1,ms).lt.0d0) jstep=0
      if(angles(2,ms).lt.0d0) kstep=0

      c1=dabs(angles(0,ms))/dx(0)
      c2=dabs(angles(1,ms))/dx(1)
      c3=dabs(angles(2,ms))/dx(2)

      do n = 1, ncells
         i=cells(0,n)+istep
         j=cells(1,n)+jstep
         k=cells(2,n)+kstep
         i_up = cells(0,n)+1-istep
         j_up = cells(1,n)+1-jstep
         k_up = cells(2,n)+1-kstep
c	print*,i,j,k,ms, i_up,j_up,k_up
         do mg=0,mg0-1
            c0=1d0/(vn(mg)*dt)+c1+c2+c3
            fi(mg,i,j,k,ms)=(qnw(mg,i,j,k,ms)
     1                      +c1*fi(mg,i_up,j,k,ms)
     2                      +c2*fi(mg,i,j_up,k,ms)
     3                      +c3*fi(mg,i,j,k_up,ms)
     4                      +qs(mg,i,j,k))/(c0+sgm_t(mg,i,j,k))
         enddo
      enddo

      return
      end   

c***********************************************************************


      subroutine setdependency(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ms0,
     &  depend0, depend1, depend2,
     &  angles)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ms0
      integer 
     &     depend0(ifirst0:ilast0+1,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2,0:ms0-1),
     &     depend1(ifirst0:ilast0,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2,0:ms0-1),
     &     depend2(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2+1,0:ms0-1)
      double precision 
     &     angles(0:3-1,0:ms0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i0, i1, i2, d0, d1, d2, d
     
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

         d2 = 0 
         if ( angles(2,d) .gt. 0d0 ) then
            d2 = 1
         else if ( angles(2,d) .lt. 0d0 ) then 
            d2 = -1
         endif

         do i2 = ifirst2, ilast2
         do i1 = ifirst1, ilast1
         do i0 = ifirst0, ilast0+1
            depend0(i0,i1,i2,d) = d0
         enddo
         enddo
         enddo

         do i2 = ifirst2, ilast2
         do i1 = ifirst1, ilast1+1
         do i0 = ifirst0, ilast0
            depend1(i0,i1,i2,d) = d1
         enddo
         enddo
         enddo

         do i2 = ifirst2, ilast2+1
         do i1 = ifirst1, ilast1
         do i0 = ifirst0, ilast0
            depend2(i0,i1,i2,d) = d2
         enddo
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************
