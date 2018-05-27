define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine preprocess_Dynamics_State(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  coord,den,cm,v,energy,vis,paq,pressure,
     &  gamma, a0, b0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1
      REAL
     &  a0,b0
c
      REAL  
     &  coord(NODE2d(ics,ice,0),0:NDIM-1),
     &  den(CELL2d(ics,ice,0)),
     &  cm(CELL2dVECG(ics,ice,icg)),
     &  v(NODE2d(ics,ice,0),0:NDIM-1),
     &  energy(CELL2d(ics,ice,0)),
     &  gamma(CELL2d(ics,ice,0)),
     &  vis(CELL2dVECG(ics,ice,icg)),
     &  paq(CELL2dVECG(ics,ice,icg)),
     &  pressure(CELL2dVECG(ics,ice,icg))
      REAL c
c
      integer i,j,k
      REAL
     & q0,x1,r1,d1,q1,q2,q3,q4,c1,s1,sukl,svkl
      REAL posvol
c
c     !计算各个网格的质量.带影像区.
      do j=ics1,ice1
      do i=ics0,ice0
         cm(i,j)=den(i,j)*posvol(
     &         coord(i,j,0),coord(i,j,1),
     &         coord(i+1,j,0),coord(i+1,j,1),
     &         coord(i+1,j+1,0),coord(i+1,j+1,1),
     &         coord(i,j+1,0),coord(i,j+1,1))
      enddo
      enddo
c
c     计算网格的压力.
      do j=ics1,ice1
      do i=ics0,ice0
         pressure(i,j) = (gamma(i,j)-1.d0)*den(i,j)*energy(i,j) 
      enddo
      enddo
c
c     计算网格的人为粘力.
c
      do j=ics1,ice1
      do i=ics0,ice0
         q1=0.0d0
         x1=coord(i+1,j+1,0)-coord(i+1,j,0)
         r1=coord(i+1,j+1,1)-coord(i+1,j,1)
         d1=sqrt(x1*x1+r1*r1)
         if(d1.ge.epsilon) then
            c1=r1/d1
            s1=x1/d1
            q1=0.5d0*( c1*(v(i+1,j+1,0)+v(i+1,j,0))
     1                -s1*(v(i+1,j+1,1)+v(i+1,j,1)))
         endif
         q2=0.0d0
         x1=coord(i+1,j+1,0)-coord(i,j+1,0)
         r1=coord(i+1,j+1,1)-coord(i,j+1,1)
         d1=sqrt(x1*x1+r1*r1)
         if(d1.ge.epsilon) then
            c1=r1/d1
            s1=x1/d1
            q2=0.5d0*(-c1*(v(i+1,j+1,0)+v(i,j+1,0))
     1                +s1*(v(i+1,j+1,1)+v(i,j+1,1)))
         endif
         q3=0.0d0
         x1=coord(i,j+1,0)-coord(i,j,0)
         r1=coord(i,j+1,1)-coord(i,j,1)
         d1=sqrt(x1*x1+r1*r1)
         if(d1.ge.epsilon) then
            c1=r1/d1
            s1=x1/d1
            q3=0.5d0*(-c1*(v(i,j+1,0)+v(i,j,0))
     1                +s1*(v(i,j+1,1)+v(i,j,1)))
         endif
         q4=0.0d0
         x1=coord(i+1,j,0)-coord(i,j,0)
         r1=coord(i+1,j,1)-coord(i,j,1)
         d1=sqrt(x1*x1+r1*r1)
         if(d1.ge.epsilon) then
            c1=r1/d1
            s1=x1/d1
            q4=0.5d0*( c1*(v(i+1,j,0)+v(i,j,0))
     1                -s1*(v(i+1,j,1)+v(i,j,1)))
         endif
c	 
         q0=q1+q2+q3+q4
c	 
         if(q0.le.0.0d0) then
           vis(i,j)=a0*den(i,j)*q0*q0

	 if(b0 > 1.0e-6) then
	   c = gamma(i,j) * pressure(i,j) / den(i,j)
           if(c < 0d0) then
c            write(*,*) i, j
c            write(*,*) "gamma = ", gamma(i,j)
c            write(*,*) "pressure = ", pressure(i,j)
c            write(*,*) "den = ", den(i,j)
c            write(*,*) "energy = ", energy(i,j)
c            write(*,*) "c < 0 when computing viscosity"
	   else 
	     c = sqrt(c)
             vis(i,j) = vis(i,j) + b0*den(i,j)*c*abs(q0)
	   endif
	 endif
	 else
	   vis(i,j)=0d0
	 endif

      enddo
      enddo
c     endif
c
c     计算网格的总压力.
      do j=ics1,ice1
      do i=ics0,ice0
         paq(i,j) = pressure(i,j) + vis(i,j)
      enddo
      enddo
c      
      return
      end 


