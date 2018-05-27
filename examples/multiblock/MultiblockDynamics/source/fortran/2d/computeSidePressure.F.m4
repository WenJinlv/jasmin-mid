define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine compute_side_pressure(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  cm,pp,
     &  side_p0, side_p1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
include(ConstantDefines.h)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
c
      REAL  
     &  pp(CELL2dVECG(ics,ice,icg)),
     &  cm(CELL2dVECG(ics,ice,icg)),
     &  side_p0(SIDE2d0VECG(ics,ice,icg)),
     &  side_p1(SIDE2d1VECG(ics,ice,icg)) 
c
      integer i,j
c
      do j=ics1-icg1,ice1+icg1
      do i=ics0,ice0+1
         side_p0(i,j) = 0.d0
         if((cm(i-1,j)+cm(i,j)).gt.1d-14)then
            side_p0(i,j)=(pp(i,j)*cm(i-1,j)+pp(i-1,j)*cm(i,j))
     &                /(cm(i-1,j)+cm(i,j))
         endif
      enddo
      enddo
c
      do j=ics1,ice1+1
      do i=ics0-icg0,ice0+icg0
         side_p1(i,j) = 0.d0
         if((cm(i,j-1)+cm(i,j)).gt.1d-14)then
            side_p1(i,j)=(pp(i,j)*cm(i,j-1)+pp(i,j-1)*cm(i,j))
     &                /(cm(i,j-1)+cm(i,j))
         endif
      enddo
      enddo
c
      return
      end 
