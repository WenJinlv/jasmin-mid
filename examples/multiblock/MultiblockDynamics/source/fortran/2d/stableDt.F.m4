define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine stable_Dt(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  coordinates,velocity,
     &  density,energy,
     &  stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dtn,rerr
      logical initial_time
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
c
      REAL  
     &  coordinates(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  velocity(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  density(CELL2dVECG(ics,ice,icg)),
     &  energy(CELL2dVECG(ics,ice,icg))
c
      integer i,j
      REAL 
     & re,dx,dr,d,c,s,er,u0,u1,u,dtc,dtkl,dt
c
c     CFL稳定性限制条件, dx < |u+-c|*dt/8.
c
c     write(*,*)"Welcome to stable.dt"
      dtkl=1.e10
      do 10 j=ics1,ice1
      do 10 i=ics0,ice0+1
         dx=coordinates(i,j+1,0)-coordinates(i,j,0)
         dr=coordinates(i,j+1,1)-coordinates(i,j,1)
         d=sqrt(dx*dx+dr*dr)
         if ( d.ge.epsilon) then 
            c=dx/d
            s=dr/d
            u0=velocity(i,  j,0)*c+velocity(i,  j,1)*s
            u1=velocity(i,j+1,0)*c+velocity(i,j+1,1)*s
            u =0.5d0*(abs(u0)+abs(u1))
            dt=0.125d0*d/(u+1.0d0)
            if(dt.le.dtkl) then
               dtkl=dt
               if(dt.le.1.d-8) then
                  write(*,"(A48,2I4,A10,E14.6,A5,E14.6)") 
     1                 "Waring : time step is very small on side-0 ",
     1                 i,j," with dt=",dtkl," dis=",d
               endif
            endif
         endif
 10   continue
      do 20 j=ics1,ice1+1
      do 20 i=ics0,ice0
         dx=coordinates(i+1,j,0)-coordinates(i,j,0)
         dr=coordinates(i+1,j,1)-coordinates(i,j,1)
         d=sqrt(dx*dx+dr*dr)
         if ( d.ge.epsilon) then
            c=dx/d
            s=dr/d
            u0=velocity(i,  j,0)*c+velocity(i,  j,1)*s
            u1=velocity(i+1,j,0)*c+velocity(i+1,j,1)*s
            u =0.5d0*(abs(u0)+abs(u1))
            dt=0.125d0*d/(u+1.0d0)
            if(dt.le.dtkl) then
               dtkl=dt
               if(dt.le.1.d-8) then
                  write(*,"(A48,2I4,A10,E14.6,A5,E14.6)") 
     1                 "Waring : time step is very small on side-1 ",
     1                 i,j," with dt=",dtkl," dis=",d
               endif
            endif
         endif
 20   continue
c
c     计算时间步长
      stabdt=dtkl
c
      return
      end 
