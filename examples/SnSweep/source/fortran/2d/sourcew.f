










c***********************************************************************
      subroutine sourceqnw(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0, ms0,
     &  dt,
     &  qnw, qw, fi,
     &  vn)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, ms0
      double precision
     &     qnw(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1),
     &     qw (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1),
     &     fi (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1)
      double precision
     &     dt,vn(0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j, ms, mg

      do ms= 0, ms0-1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg = 0, mg0-1
         qnw(mg,i,j,ms)=fi(mg,i,j,ms)/(vn(mg)*dt)
     &                   +qw(mg,i,j,ms)
      enddo
      enddo
      enddo
      enddo

      return
      end

c***********************************************************************

      subroutine sourcew(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0, ms0,
     &  qw,
     &  dx, xleng,
     &  angles,
     &  vn,
     &  ysgm_s,
     &  ysgm_t,
     &  yd,
     &  am, bm)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, ms0
      double precision
     &     qw(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1)
      double precision
     &     dx(0:2-1), xleng(0:2-1),
     &     angles(0:2-1,0:ms0-1),
     &     vn(0:mg0-1),ysgm_s(0:mg0-1),ysgm_t(0:mg0-1)
      double precision 
     &     yd,
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j, ms,mg
      integer l0,l1
      double precision    skes, smiu
      double precision    aa, bb
     
      
      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo


      do ms= 0, ms0-1
         l0 = 0
         l1 = 0
         skes = angles(0,ms)
         if ( skes .lt. 0.d0 ) then
            skes = -skes
            l0 = 1
         endif
         smiu = angles(1,ms)
         if ( smiu .lt. 0.d0 ) then
            smiu = -smiu
            l1 = 1
         endif
         do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0+1
            aa = am(i,l0)*bm(j,l1)
            bb = skes*bm(j,l1)/xleng(0)
     &          +smiu*am(i,l0)/xleng(1)
            do mg = 0, mg0-1
               qw(mg,i,j,ms) = yd*( - ysgm_s(mg)
     &                              + ysgm_t(mg)*aa + vn(mg)*bb )
            enddo
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************

      subroutine diffsolution(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0, ms0,
     &  fi,
     &  dx, xleng,
     &  yd, vn, angles,
     &  am, bm)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, ms0
      double precision
     &     fi(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1)
      double precision 
     &     dx(0:2-1), xleng(0:2-1),
     &     yd, vn(0:mg0-1),
     &     angles(0:2-1,0:ms0-1),
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j, ms,mg, l0,l1
      double precision aa, solution, err

      
      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo


      do ms= 0, ms0-1
         l0 = 0
         l1 = 0
         if ( angles(0,ms).lt. 0.d0 ) l0 = 1
         if ( angles(1,ms).lt. 0.d0 ) l1 = 1
         do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0+1
            aa = yd*am(i,l0)*bm(j,l1)
            do mg = 0, mg0-1
               solution = vn(mg)*aa
	       if ( solution .le. 1.d0 ) then
                   err = abs(fi(mg,i,j,ms)-solution)
               else 
                   err = abs((fi(mg,i,j,ms)-solution)/solution)
               endif
               if ( err.gt.1.d-4) then
                  print*, "ERROR: ", 
     &                    i,j,mg,ms,solution,fi(mg,i,j,ms)
                  stop
               endif
            enddo
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************
