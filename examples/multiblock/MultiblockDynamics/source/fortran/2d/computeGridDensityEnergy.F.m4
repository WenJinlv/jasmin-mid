define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine compute_GridDensityEnergy(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  coord,coord_new,
     &  density, density_new,
     &  energy, energy_new,
     &  pressure,
     &  paq,vel,viscosity,
     &  gamma,
     &  dt)
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
     &  coord(NODE2d(ics,ice,0),0:NDIM-1),
     &  coord_new(NODE2d(ics,ice,0),0:NDIM-1),
     &  density(CELL2d(ics,ice,0)),
     &  density_new(CELL2d(ics,ice,0)),
     &  energy(CELL2d(ics,ice,0)),
     &  energy_new(CELL2d(ics,ice,0)),
     &  gamma(CELL2d(ics,ice,0)),
     &  vel(NODE2d(ics,ice,0),0:NDIM-1),
     &  pressure(CELL2dVECG(ics,ice,icg)),
     &  paq(CELL2dVECG(ics,ice,icg)),
     &  viscosity(CELL2dVECG(ics,ice,icg))
c
      integer i,j
      REAL 
     &  area, posvol,
     &  v, vn,
     &  de_half, dDen_half, dP_half
c
c
      do 30 j=ics1,ice1+1
      do 30 i=ics0,ice0+1

c 更新网格坐标
         coord_new(i,j,0)=coord(i,j,0)+dt*vel(i,j,0)
         coord_new(i,j,1)=coord(i,j,1)+dt*vel(i,j,1)
         if(coord_new(i,j,1).le.0.0d0) then
            coord_new(i,j,1)=0.0d0
            vel(i,j,1)=0.0d0
         endif
 30   continue


c 更新网格密度
      do 40 j=ics1,ice1
      do 40 i=ics0,ice0

         area= 0.5d0*( (coord_new(i+1,j+1,0)-coord_new(i,j,0))
     &                  *(coord_new(i,j+1,1)-coord_new(i+1,j,1))
     &                  +(coord_new(i+1,j,0)-coord_new(i,j+1,0))
     &                  *(coord_new(i+1,j+1,1)-coord_new(i,j,1)) )

         if(area.lt.0.0d0) then
            print *,"k,l,negative area=",i,j,area,
     &     "-------------------------------------------------------",
     &      coord(i,j,0),coord(i,j,1),
     &      coord(i+1,j,0),coord(i+1,j,1),
     &      coord(i+1,j+1,0),coord(i+1,j+1,1),
     &      coord(i,j+1,0),coord(i,j+1,1),
     &    "******************************************************",
     &      paq(i,j),viscosity(i,j),dt,
     &      vel(i,j,0),vel(i,j,1),
     &      vel(i+1,j,0),vel(i+1,j,1),
     &      vel(i+1,j+1,0),vel(i+1,j+1,1),
     &      vel(i,j+1,0),vel(i,j+1,1),
     &      "=======================================================",
     &      coord_new(i,j,0),coord_new(i,j,1),
     &      coord_new(i+1,j,0),coord_new(i+1,j,1),
     &      coord_new(i+1,j+1,0),coord_new(i+1,j+1,1),
     &      coord_new(i,j+1,0),coord_new(i,j+1,1)
         endif

         vn=posvol(coord(i,j,0),coord(i,j,1),
     1             coord(i+1,j,0),coord(i+1,j,1),
     1             coord(i+1,j+1,0),coord(i+1,j+1,1),
     1             coord(i,j+1,0),coord(i,j+1,1))
         v =posvol(coord_new(i,j,0),coord_new(i,j,1),
     2             coord_new(i+1,j,0),coord_new(i+1,j,1),
     2             coord_new(i+1,j+1,0),coord_new(i+1,j+1,1),
     2             coord_new(i,j+1,0),coord_new(i,j+1,1))

         if(v .le. 0.0d0) then
             write(*,"('negative volume of mesh(',I4,',',I4,')')")
     &                i, j
             write(*,"('vol(i,j) = ', E16.5)") v
         endif
         density_new(i,j)= density(i,j)*vn/v
   40 continue
    
c 更新网格能量

c//能量方程
c//两步格式计算内能
c//参见李德元《二维非定常流体力学数值方法》p139
      do 50 j=ics1,ice1
      do 50 i=ics0,ice0

c        //第一步：预估
         de_half = energy(i,j) - 0.5d0 * paq(i,j)
     &             *(1.d0/density_new(i,j) - 1.d0/density(i,j))

c        //第二步：校正
         dDen_half = (density_new(i,j) + density(i,j)) * 0.5d0
         dP_half   = (gamma(i,j)-1.d0)*dDen_half*de_half
     &                          + viscosity(i,j) 

         energy_new(i,j) = energy(i,j) - dP_half *
     &               (1.d0 / density_new(i,j) - 1.d0 / density(i,j))

        if (energy_new(i,j) < 0.d0)then
           write(*,*) "energy_new < 0 on cell",i,"   ",j
        end if 

   50 continue
c
       return
      end 
