define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    求累加和.
c
c***********************************************************************
      subroutine advancing(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,dt,
     &  nc, source, uo, un)
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
      REAL dt,
     &     nc(NODE2d(ifirst,ilast,0)),
     &     source(CELL2d(ifirst,ilast,0)),
     &     uo(CELL2dVECG(ifirst,ilast,gcw)),
     &     un(CELL2d(ifirst,ilast,0))
c
      REAL dx, cfl, kapp
      integer J,K
      logical is_print
c
c***********************************************************************     
c

      is_print = .false.
      cfl = 2.0d-1
      kapp = 1.0d0

      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0
         un(j,k)=   (1.0-4.0*cfl*kapp)*uo(j,k)+
     &             kapp*cfl*(uo(j,k-1)+uo(j,k+1)
     &             +uo(j-1,k)+uo(j+1,k))
     &             +dt*source(j,k)
c         nc(j  ,k)  =   nc(j  ,k)  +uo(j,k)
c         nc(j+1,k)  =   nc(j+1,k)  +uo(j,k)
c         nc(j,k+1)  =   nc(j,k+1)  +uo(j,k)
c         nc(j+1,k+1)=   nc(j+1,k+1)+uo(j,k)
      ENDDO
      ENDDO
c
      return
      end


