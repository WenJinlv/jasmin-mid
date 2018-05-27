define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    求累加和.
c
c***********************************************************************
      subroutine filling_phybdr(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  ifirst_f0,ilast_f0,ifirst_f1,ilast_f1,
     &  var_current, var_scratch)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ifirst_f0,ilast_f0,ifirst_f1,ilast_f1
      integer gcw0,gcw1
c
c variables in 2d indexed
c
      REAL 
     &     var_scratch(CELL2dVECG(ifirst,ilast,gcw)),
     &     var_current(CELL2d(ifirst,ilast,0))
c
      integer i,j, size_f_0, size_f_1
      logical is_print
c
c***********************************************************************     
c

      size_f_0 = ilast_f0-ifirst_f0+1
      size_f_1 = ilast_f1-ifirst_f1+1

      if (size_f_0 .eq. 1) then

         i = ifirst_f0
c        left:
         if (i+1.eq.ifirst0) then
            DO j=ifirst_f1,ilast_f1
               var_scratch(i,j) = var_current(ifirst0,j) 
            ENDDO

c        right:
         else if (i-1.eq.ilast0) then
            DO j=ifirst_f1,ilast_f1
               var_scratch(i,j) = var_current(ilast0,j) 
            ENDDO

         else 
            print*, "left or right phybdry error."
            stop
         endif

      else if (size_f_1 .eq. 1) then

         j = ifirst_f1
c        bottom:
         if (j+1.eq.ifirst1) then
            DO i=ifirst_f0,ilast_f0
               var_scratch(i,j) = var_current(i,ifirst1) 
            ENDDO

c        up:
         else if (j-1.eq.ilast1) then
            DO i=ifirst_f0,ilast_f0
               var_scratch(i,j) = var_current(i,ilast1) 
            ENDDO

         else 
            print*, "bottom or up phybdry error."
            stop
         endif


      else 
         print*, "phybdry error."
         stop
      endif
c
      return
      end


