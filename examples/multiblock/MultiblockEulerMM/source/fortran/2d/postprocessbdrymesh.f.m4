define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    处理边界处的网格.
c
c***********************************************************************
      subroutine postprocessbdrymesh(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  igfirst0,iglast0,igfirst1,iglast1,
     &  type,loc,xo,yo,xn,yn)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,
     &        igfirst0,iglast0,igfirst1,iglast1
      integer gcw0,gcw1,type,loc
c
c variables in 2d indexed
c
      REAL
     &     xo(NODE2dVECG(ifirst,ilast,gcw)),
     &     yo(NODE2dVECG(ifirst,ilast,gcw)),
     &     xn(NODE2d(ifirst,ilast,0)),
     &     yn(NODE2d(ifirst,ilast,0))
c
      integer i,j,l
c
c***********************************************************************     
C
      if(type.eq.1) then
         if(loc.eq.0) then
            do j=max0(igfirst1,ifirst1),min0(iglast1+1,ilast1+1)
               xn(iglast0+1,j) =  xo(iglast0+1,j)
               yn(iglast0+1,j) =  yo(iglast0+1,j)
     &                           +(yn(iglast0+2,j)-yo(iglast0+2,j))
            enddo
         else if(loc.eq.1) then
            do j=max0(igfirst1,ifirst1),min0(iglast1+1,ilast1+1)
               xn(igfirst0,j) =  xo(igfirst0,j)
               yn(igfirst0,j) =  yo(igfirst0,j)
     &                           +(yn(igfirst0-1,j)-yo(igfirst0-1,j))
            enddo
         else if(loc.eq.2) then
            do i=max0(igfirst0,ifirst0),min0(iglast0+1,ilast0+1)
               xn(i,iglast1+1) =  xo(i,iglast1+1)
     &                           +(xn(i,iglast1+2)-xo(i,iglast1+2))
               yn(i,iglast1+1) =  yo(i,iglast1+1)
            enddo
         else 
            do i=max0(igfirst0,ifirst0),min0(iglast0+1,ilast0+1)
               xn(i,igfirst1) =  xo(i,igfirst1)
     &                           +(xn(i,igfirst1-1)-xo(i,igfirst1-1))
               yn(i,igfirst1) =  yo(i,igfirst1)
            enddo
         endif
      endif
C
      return
      end
