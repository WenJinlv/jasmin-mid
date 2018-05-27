

c***********************************************************************
      subroutine sourceqnw(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0, ms0,
     &  dt,
     &  qnw, qw, fi,
     &  vn)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      double precision
     &     qnw(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1),
     &     qw (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1),
     &     fi (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1)
      double precision
     &     dt,vn(0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j,k, ms, mg

      do ms= 0, ms0-1
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg = 0, mg0-1
         qnw(mg,i,j,k,ms)=fi(mg,i,j,k,ms)/(vn(mg)*dt)
     &                   +qw(mg,i,j,k,ms)
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

c***********************************************************************

      subroutine sourcew(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0, ms0,
     &  qw,
     &  dx, xleng,
     &  angles,
     &  vn,
     &  ysgm_s,
     &  ysgm_t,
     &  yd,
     &  am, bm, cm)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      double precision
     &     qw(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1)
      double precision
     &     dx(0:3-1), xleng(0:3-1),
     &     angles(0:3-1,0:ms0-1),
     &     vn(0:mg0-1),ysgm_s(0:mg0-1),ysgm_t(0:mg0-1)
      double precision 
     &     yd,
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1),
     &     cm(ifirst2:ilast2+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j,k, ms,mg
      integer l0,l1,l2
      double precision    skes, smiu, seta
      double precision    aa, bb
     
      
      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo
      do k = ifirst2, ilast2+1
         cm(k,0) = k*dx(2)/xleng(2)
         cm(k,1) = 1d0 - cm(k,0)
      enddo


      do ms= 0, ms0-1
         l0 = 0
         l1 = 0
         l2 = 0
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
         seta = angles(2,ms)
         if ( seta .lt. 0.d0 ) then
            seta = -seta
            l2 = 1
         endif
         do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0+1
            aa = am(i,l0)*bm(j,l1)*cm(k,l2)
            bb = skes*bm(j,l1)*cm(k,l2)/xleng(0)
     &          +smiu*cm(k,l2)*am(i,l0)/xleng(1)
     &          +seta*am(i,l0)*bm(j,l1)/xleng(2)
            do mg = 0, mg0-1
               qw(mg,i,j,k,ms) = yd*( - ysgm_s(mg)
     &                                + ysgm_t(mg)*aa + vn(mg)*bb )
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************

      subroutine diffsolution(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0, ms0,
     &  fi,
     &  dx, xleng,
     &  yd, vn, angles,
     &  am, bm, cm)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      double precision
     &     fi(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:ms0-1)
      double precision 
     &     dx(0:3-1), xleng(0:3-1),
     &     yd, vn(0:mg0-1),
     &     angles(0:3-1,0:ms0-1),
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1),
     &     cm(ifirst2:ilast2+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j,k, ms,mg, l0,l1,l2
      double precision aa, err

      
      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo
      do k = ifirst2, ilast2+1
         cm(k,0) = k*dx(2)/xleng(2)
         cm(k,1) = 1d0 - cm(k,0)
      enddo


      do ms= 0, ms0-1
         l0 = 0
         l1 = 0
         l2 = 0
         if ( angles(0,ms).lt. 0.d0 ) l0 = 1
         if ( angles(1,ms).lt. 0.d0 ) l1 = 1
         if ( angles(2,ms).lt. 0.d0 ) l2 = 1
         do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0+1
            aa = yd*am(i,l0)*bm(j,l1)*cm(k,l2)
            do mg = 0, mg0-1
               err = 0d0
               if ( fi(mg,i,j,k,ms) .ne. 0d0 ) then
                  err = (fi(mg,i,j,k,ms)-vn(mg)*aa)/fi(mg,i,j,k,ms)
               endif
               if ( err .gt. 1.d-5) then
                  print*, "ERROR: ", 
     &                    i,j,k,mg,ms,fi(mg,i,j,k,ms),vn(mg)*aa,err
                  stop
               endif
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end   

c***********************************************************************
