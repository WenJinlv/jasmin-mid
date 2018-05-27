define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine compute_Velocity(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  vel,vel_new,sp,dt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
      REAL
     &  dt
c
      REAL
     &  vel(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  vel_new(NODE2d(ics,ice,0),0:NDIM-1),
     &  sp(NODE2dVECG(ics,ice,icg),0:NDIM-1)
c
      integer i,j
c
c     更新网格片内部结点的速度.
c

      do j=ics1,ice1+1
      do i=ics0,ice0+1
         vel_new(i,j,0)=vel(i,j,0)-dt*sp(i,j,0)
         vel_new(i,j,1)=vel(i,j,1)-dt*sp(i,j,1)
         if(abs(vel_new(i,j,0)).lt.1.d-6) vel_new(i,j,0) = 0.0d0
         if(abs(vel_new(i,j,1)).lt.1.d-6) vel_new(i,j,1) = 0.0d0
c       write(*,*)i,j,vel(i,j,0),vel(i,j,1)
c       write(*,*)i,j,sp(i,j,0),sp(i,j,1)
c       write(*,*)i,j,vel_new(i,j,0),vel_new(i,j,1)
      enddo
      enddo
c
      return
      end 
