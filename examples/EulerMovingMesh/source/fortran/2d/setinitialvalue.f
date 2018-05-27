










c***********************************************************************
c
c    赋初值.
c
c***********************************************************************
      subroutine setinitialvalue(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,domain_id,gamma,
     &  uval)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1,domain_id
      double precision
     &     gamma
c variables in 2d cell indexed
      double precision
     &     uval(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,2+2)
c
      integer i,j
      double precision
     &  rho, vel0, vel1, eng
c
c***********************************************************************     
c
      if (domain_id==0) then 
          rho  = 0.8d0
          vel0 = 0.1d0
          vel1 = 0.1d0
          eng  = 0.4d0
      else if(domain_id==1) then
          rho  = 0.5197d0
          vel0 = 0.1d0
          vel1 = -0.6259d0
          eng  = 0.4d0
      else if(domain_id==2) then
          rho  = 0.5197d0
          vel0 = -0.6259d0
          vel1 = 0.1d0
          eng  = 0.4d0
      else if(domain_id==3) then
          rho  = 1.0d0
          vel0 = 0.1d0
          vel1 = 0.1d0
          eng  = 1.0d0
      else
          stop "wrong domain_id : "
      endif
c
      do j=ifirst1, ilast1
      do i=ifirst0, ilast0
         uval(i,j,1) = rho
         uval(i,j,2) = uval(i,j,1)*vel0
         uval(i,j,3) = uval(i,j,1)*vel1
         uval(i,j,4) = eng/(gamma-1.0)
     &                + 0.5*(uval(i,j,2)**2+uval(i,j,3)**2)
     &                  /uval(i,j,1)
      enddo
      enddo
c
      return
      end
