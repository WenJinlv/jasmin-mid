define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine postprocess_Velocity(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  bns0,bne0,bns1,bne1,
     &  coord,vel,
     &  location,type)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
include(ConstantDefines.h)dnl
include(BoundaryDefines.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  bns0,bne0,bns1,bne1,
     &  location,type,dir
c
      REAL  
     &  coord(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  vel(NODE2d(ics,ice,0),0:NDIM-1)
c
      integer k,l,kst,kend,lst,lend,ka1,la1,km1,lm1
      REAL    
     &  dx,dr,dd,cosa,sina,x1,r1,x2,r2,x3,r3,x4,r4,u1,v1,
     &  u2,v2,u3,v3,au,av,dx1,dx2,dr1,dr2,dx3,dr3,
     &  a2,b2,c2,a1,b1,a,al,value,
     &  x10,r10,x20,r20,u10,v10,u20,v20,
     &  fm,fx,x5,fr,r5,d10,d20,u5,v5
     &  ,vnew0,vnew1
c
c     更新网格结点的速度.
c
      kst =bns0
      kend=bne0
      lst =bns1
      lend=bne1
      if(location.eq.XLO) then
         kst  = ics0
         kend = kst
         lst  = ics1
         lend = ice1+1
         dir = 0
      else if(location.eq.XHI) then
         kst  = ice0+1
         kend = kst
         lst  = ics1
         lend = ice1+1
         dir = 0
      else if(location.eq.YLO) then
         kst  = ics0
         kend = ice0+1
         lst  = ics1
         lend = lst
         dir = 1
      else if(location.eq.YHI) then
         kst  = ics0
         kend = ice0+1
         lst  = ice1+1
         lend = lst
         dir = 1
      endif
c
      if(type.eq.SYMM_BDY2D_COND_TYPE .or.
     &   type.eq.WALL_BDY2D_COND_TYPE) then
c
         do 30 l=lst,lend
         do 30 k=kst,kend
            ka1=k+dir
            la1=l+1-dir
            dx=coord(ka1,la1,0)-coord(k,l,0)
            dr=coord(ka1,la1,1)-coord(k,l,1)
            dd=sqrt(dx*dx+dr*dr)
            cosa=dx/dd
            sina=dr/dd
            if (abs(cosa).lt.epsilon) then 
                cosa=0d0
                sina=1d0
            endif
            if (abs(sina).lt.epsilon) then 
                cosa=1d0
                sina=0d0
            endif
            vnew0=(vel(k,l,0)*cosa+vel(k,l,1)*sina)*cosa
            vnew1=(vel(k,l,0)*cosa+vel(k,l,1)*sina)*sina
            vel(k,l,0) = vnew0
            vel(k,l,1) = vnew1
            if(abs(vel(k,l,0)).lt.1.d-6) vel(k,l,0) = 0.0d0
            if(abs(vel(k,l,1)).lt.1.d-6) vel(k,l,1) = 0.0d0
 30      continue
       endif

       if(type.eq.ZERO_BDY2D_COND_TYPE) then
         do l=lst,lend
         do k=kst,kend
            vel(k,l,0) = 0.0d0
            vel(k,l,1) = 0.0d0
	 enddo
	 enddo
      endif

       return
       end 
