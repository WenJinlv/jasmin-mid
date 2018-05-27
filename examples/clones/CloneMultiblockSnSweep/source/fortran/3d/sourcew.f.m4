define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim3d.i)dnl
include(m4func.i)dnl
c***********************************************************************
      subroutine sourceqnw(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0, ms0,
     &  dt,
     &  qnw, qw, fi,
     &  vn,
     &  mg_begin, mg_end, mg_end2, loc2glob)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2
      integer loc2glob(mg_begin: mg_end2),ig

      REAL
     &     qnw(mg_begin:mg_end2,NODE3d(ifirst,ilast,0),0:ms0-1),
     &     qw (mg_begin:mg_end2,NODE3d(ifirst,ilast,0),0:ms0-1),
     &     fi (mg_begin:mg_end2,NODE3d(ifirst,ilast,0),0:ms0-1)
      REAL
     &     dt,vn(0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j,k, ms, mg

      do ms= 0, ms0-1
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
      do mg = mg_begin, mg_end
         ig = loc2glob( mg )
         qnw(mg,i,j,k,ms)=fi(mg,i,j,k,ms)/(vn(ig)*dt)
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
     &  am, bm, cm,
     &  mg_begin, mg_end, mg_end2, loc2glob)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2
      integer loc2glob(mg_begin:mg_end2)
      integer ig

      REAL
     &  qw(mg_begin:mg_end2, NODE3d(ifirst,ilast,0),0:ms0-1)
      REAL
     &     dx(0:NDIM-1), xleng(0:NDIM-1),
     &     angles(0:NDIM-1,0:ms0-1),
     &     vn(0:mg0-1),ysgm_s(0:mg0-1),ysgm_t(0:mg0-1)
      REAL 
     &     yd,
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1),
     &     cm(ifirst2:ilast2+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j,k, ms,mg
      integer l0,l1,l2
      REAL    skes, smiu, seta
      REAL    aa, bb
     
      CAL_AM_BM_CM

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
            do mg = mg_begin, mg_end
               ig = loc2glob( mg )
               qw(mg,i,j,k,ms) = yd*( - ysgm_s(ig)
     &                + ysgm_t(ig)*aa + vn(ig)*bb )
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
     &  am, bm, cm,
     &  mg_begin, mg_end, mg_end2, loc2glob)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0, ms0
      integer mg_begin, mg_end, mg_end2
      integer loc2glob(mg_begin:mg_end2)

      REAL
     &     fi(mg_begin:mg_end, NODE3d(ifirst,ilast,0),0:ms0-1)
      REAL 
     &     dx(0:NDIM-1), xleng(0:NDIM-1),
     &     yd, vn(0:mg0-1),
     &     angles(0:NDIM-1,0:ms0-1),
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1),
     &     cm(ifirst2:ilast2+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      integer i,j,k, ms,mg, l0,l1,l2, ig
      REAL aa, err

      CAL_AM_BM_CM

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
            do mg = mg_begin, mg_end
               ig = loc2glob(mg)
               err = 0d0
               if ( fi(mg,i,j,k,ms) .ne. 0d0 ) then
                  err = (fi(mg,i,j,k,ms)-vn(ig)*aa)/fi(mg,i,j,k,ms)
               endif
               if ( err .gt. 1.d-5) then
                  print*, "ERROR: ", 
     &                    i,j,k,mg,ms,fi(mg,i,j,k,ms),vn(ig)*aa,err
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
