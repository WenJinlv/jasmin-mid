define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine initialize_velocity(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  coord,vel_cur,vv)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
c
      REAL
     &  vel_cur(NODE2d(ics,ice,0),0:NDIM-1),
     &  coord(NODE2d(ics,ice,0),0:NDIM-1),
     &  vv
c
      integer i,j
      REAL x,y,d
c
c     更新网格片内部结点的速度.
c
      do j=ics1,ice1+1
      do i=ics0,ice0+1
         x = coord(i,j,0)
         y = coord(i,j,1)
         if(abs(x) < 1.0e-8 .and. abs(y) < 1.0e-8)  then
           vel_cur(i,j,0)=0d0
           vel_cur(i,j,1)=0d0
	 elseif(abs(x) < 1.0e-8) then
           vel_cur(i,j,0)=0d0
           vel_cur(i,j,1)=-vv
	 elseif(abs(y) < 1.0e-8) then
           vel_cur(i,j,0)=-vv
           vel_cur(i,j,1)=0d0
	 else
	   d=sqrt(x*x+y*y)
           vel_cur(i,j,0)=-vv * (x/d)
           vel_cur(i,j,1)=-vv * (y/d)
         endif
      enddo
      enddo
c
      return
      end 


      subroutine zeroleft_velocity(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  vel_cur)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
c
      REAL
     &  vel_cur(NODE2d(ics,ice,0),0:NDIM-1)
c
      integer i,j
c
      do j=ics1,ice1+1
      do i=ics0,ics0
           vel_cur(i,j,0)=0d0
           vel_cur(i,j,1)=0d0
      enddo
      enddo
c
      return
      end 


      subroutine zerobottom_velocity(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  vel_cur)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
c
      REAL
     &  vel_cur(NODE2d(ics,ice,0),0:NDIM-1)
c
      integer i,j
c
      do j=ics1,ics1
      do i=ics0,ice0+1
           vel_cur(i,j,0)=0d0
           vel_cur(i,j,1)=0d0
      enddo
      enddo
c

      return
      end 


