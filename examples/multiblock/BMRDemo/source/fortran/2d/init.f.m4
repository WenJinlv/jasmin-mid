define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    求累加和.
c
c***********************************************************************
      subroutine init(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  source)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
c
c variables in 2d indexed
c
      REAL
     &     source(CELL2d(ifirst,ilast,0))
c
      integer J,K
      logical is_print
c
c***********************************************************************     
c

      is_print = .false.

      k = ifirst1
c      DO k=ifirst1,ilast1
      DO J=ifirst0,ilast0
         source(j,k)=1.0d4
c      ENDDO
      ENDDO

c
      return
      end


