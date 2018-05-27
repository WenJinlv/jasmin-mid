define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine alternal_node_coupler(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  ibs0,ibe0,ibs1,ibe1,
     &  vel,vel_tmp,xi)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
include(BoundaryDefines.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  ibs0,ibe0,ibs1,ibe1
      REAL
     &  xi
c
      REAL
     &  vel(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  vel_tmp(NODE2dVECG(ics,ice,icg),0:NDIM-1)
c
      integer i,j
      REAL eta

c
c     修正网格片所有结点的速度.
      eta = 0.5d0*(1.0d0+xi)
      do j=ibs1,ibe1
      do i=ibs0,ibe0
         vel_tmp(i,j,0) = vel_tmp(i,j,0)
     &                      +eta*(vel(i,j+1,0)+vel(i+1,j,0))
     &                      -xi*vel(i+1,j+1,0)-vel(i,j,0) 
         vel_tmp(i,j,1) = vel_tmp(i,j,1)
     &                      +eta*(vel(i,j+1,1)+vel(i+1,j,1))
     &                      -xi*vel(i+1,j+1,1)-vel(i,j,1)

         vel_tmp(i+1,j,0) = vel_tmp(i+1,j,0)
     &                       +eta*(vel(i,j,0)+vel(i+1,j+1,0))
     &                       -xi*vel(i,j+1,0)-vel(i+1,j,0) 
         vel_tmp(i+1,j,1) = vel_tmp(i+1,j,1)
     &                       +eta*(vel(i,j,1)+vel(i+1,j+1,1))
     &                       -xi*vel(i,j+1,1)-vel(i+1,j,1) 

         vel_tmp(i+1,j+1,0) = vel_tmp(i+1,j+1,0)
     &                      +eta*(vel(i,j+1,0)+vel(i+1,j,0))
     &                      -xi*vel(i,j,0)-vel(i+1,j+1,0) 
         vel_tmp(i+1,j+1,1) = vel_tmp(i+1,j+1,1)
     &                      +eta*(vel(i,j+1,1)+vel(i+1,j,1))
     &                      -xi*vel(i,j,1)-vel(i+1,j+1,1)

         vel_tmp(i,j+1,0) = vel_tmp(i,j+1,0)
     &                       +eta*(vel(i,j,0)+vel(i+1,j+1,0))
     &                       -xi*vel(i+1,j,0)-vel(i,j+1,0) 
         vel_tmp(i,j+1,1) = vel_tmp(i,j+1,1)
     &                       +eta*(vel(i,j,1)+vel(i+1,j+1,1))
     &                       -xi*vel(i+1,j,1)-vel(i,j+1,1) 
      enddo
      enddo
c
      return
      end 
