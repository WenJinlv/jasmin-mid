define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
      subroutine compute_Speedup_IGA(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  ibs0,ibe0,ibs1,ibe1,
     &  coord,den,pp,
     &  side_p0,side_p1,
     &  sp)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
include(BoundaryDefines.i)dnl
include(ConstantDefines.h)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  ibs0,ibe0,ibs1,ibe1
      integer
     & is,ie,js,je
c
      REAL  
     &  coord(NODE2dVECG(ics,ice,icg),0:NDIM-1),
     &  den(CELL2dVECG(ics,ice,icg)),
     &  pp(CELL2dVECG(ics,ice,icg)),
     &  side_p0(SIDE2d0VECG(ics,ice,icg)),
     &  side_p1(SIDE2d1VECG(ics,ice,icg)),
     &  sp(NODE2dVECG(ics,ice,icg),0:NDIM-1)
c
      integer i,j,k
      REAL  dden,pres,posvol,negvol
      REAL  x(0:2*NDIM),r(0:2*NDIM),vol0(2*NDIM),
     &          x0(0:2*NDIM),r0(0:2*NDIM),
     &          p0(2*NDIM),ps(4,2*NDIM),
     &          sx(2*NDIM),sr(2*NDIM),omg(2*NDIM)

c
c        !以网格单元为索引,计算各个单元对其四个结点加速度的贡献,并累加到结点上.
c

      do j=ibs1,ibe1
      do i=ibs0,ibe0
        dden=den(i,j)
        pres=pp(i,j)

        !临时存储四个网格结点和中心点的坐标
        x(1)=coord(i , j,0)
        x(2)=coord(i+1,j,0)
        x(3)=coord(i+1,j+1,0)
        x(4)=coord(i  ,j+1,0)
        x(0)=0.25d0*(x(1)+x(2)+x(3)+x(4))
        r(1)=coord(i , j,1) 
        r(2)=coord(i+1,j,1) 
        r(3)=coord(i+1,j+1 ,1) 
        r(4)=coord(i  ,j+1 ,1)
        r(0)=0.25d0*(r(1)+r(2)+r(3)+r(4))
c
        !分解网格单元为四个子单元,分别对应四个结点,并为子单元和网格单元计算体积.
        x0(1)=0.5d0*(x(1)+x(2))
        x0(2)=0.5d0*(x(2)+x(3))
        x0(3)=0.5d0*(x(3)+x(4))
        x0(4)=0.5d0*(x(4)+x(1))
        x0(0)=0.25d0*(x(1)+x(2)+x(3)+x(4))
        r0(1)=0.5d0*(r(1)+r(2))
        r0(2)=0.5d0*(r(2)+r(3))
        r0(3)=0.5d0*(r(3)+r(4))
        r0(4)=0.5d0*(r(4)+r(1))
        r0(0)=0.25d0*(r(1)+r(2)+r(3)+r(4))
        vol0(1)=negvol(x(1),r(1),x0(1),r0(1),x(0),r(0),x0(4),r0(4))
        vol0(2)=negvol(x0(1),r0(1),x(2),r(2),x0(2),r0(2),x(0),r(0))
        vol0(3)=negvol(x(0),r(0),x0(2),r0(2),x(3),r(3),x0(3),r0(3))
        vol0(4)=negvol(x0(4),r0(4),x(0),r(0),x0(3),r0(3),x(4),r(4))
c
            !为每个结点,计算子单元四条边上的压力.每个子单元的结点,以左下角结点为起点.
            ! 第1个结点的子单元
            ps(1,1)=side_p1(i,j)
            ps(2,1)=pp(i,j)
            ps(3,1)=pp(i,j)
            ps(4,1)=side_p0(i,j)
            ! 第2个结点的子单元
            ps(1,2)=ps(1,1)
            ps(2,2)=side_p0(i+1,j)
            ps(3,2)=pp(i,j)
            ps(4,2)=pp(i,j)
            !第3个结点的子单元
            ps(1,3)=pp(i,j)
            ps(2,3)=ps(2,2)
            ps(3,3)=side_p1(i,j+1)
            ps(4,3)=pp(i,j)
            !第4个结点的子单元
            ps(1,4)=pp(i,j)
            ps(2,4)=pp(i,j)
            ps(3,4)=ps(3,3)
            ps(4,4)=ps(4,1)      
c
c           计算四个结点上的加速度.
            !X方向的压力分量
            sx(1)= ps(1,1)*(r0(1)-r(1)) *(r0(1)+r(1))
     &            +ps(2,1)*(r0(0)-r0(1))*(r0(0)+r0(1))
     &            +ps(3,1)*(r0(4)-r0(0))*(r0(4)+r0(0))
     &            +ps(4,1)*(r(1)-r0(4)) *(r(1)+r0(4))
            sx(2)= ps(1,2)*(r(2)-r0(1)) *(r(2)+r0(1))
     &            +ps(2,2)*(r0(2)-r(2)) *(r0(2)+r(2))
     &            +ps(3,2)*(r0(0)-r0(2))*(r0(0)+r0(2))
     &            +ps(4,2)*(r0(1)-r0(0))*(r0(1)+r0(0))
            sx(3)= ps(1,3)*(r0(2)-r0(0))*(r0(2)+r0(0))
     &            +ps(2,3)*(r(3)-r0(2)) *(r(3)+r0(2))
     &            +ps(3,3)*(r0(3)-r(3)) *(r0(3)+r(3))
     &            +ps(4,3)*(r0(0)-r0(3))*(r0(0)+r0(3))
            sx(4)= ps(1,4)*(r0(0)-r0(4))*(r0(0)+r0(4))
     &            +ps(2,4)*(r0(3)-r0(0))*(r0(3)+r0(0))
     &            +ps(3,4)*(r(4)-r0(3)) *(r(4)+r0(3))
     &            +ps(4,4)*(r0(4)-r(4)) *(r0(4)+r(4))
            !R方向的压力分量
          sr(1)=  -ps(1,1)*(x0(1)-x(1)) *(r0(1)+r(1))
     &            -ps(2,1)*(x0(0)-x0(1))*(r0(0)+r0(1))
     &            -ps(3,1)*(x0(4)-x0(0))*(r0(4)+r0(0))
     &            -ps(4,1)*(x(1)-x0(4)) *(r(1)+r0(4))
          omg(1)= -(x0(1)-x(1)) *(r0(1)+r(1))
     &            -(x0(0)-x0(1))*(r0(0)+r0(1))
     &            -(x0(4)-x0(0))*(r0(4)+r0(0))
     &            -(x(1)-x0(4)) *(r(1)+r0(4))
          sr(2)=  -ps(1,2)*(x(2)-x0(1)) *(r(2)+r0(1))
     &            -ps(2,2)*(x0(2)-x(2)) *(r0(2)+r(2))
     &            -ps(3,2)*(x0(0)-x0(2))*(r0(0)+r0(2))
     &            -ps(4,2)*(x0(1)-x0(0))*(r0(1)+r0(0))
          omg(2)= -(x(2)-x0(1)) *(r(2)+r0(1))
     &            -(x0(2)-x(2)) *(r0(2)+r(2))
     &            -(x0(0)-x0(2))*(r0(0)+r0(2))
     &            -(x0(1)-x0(0))*(r0(1)+r0(0))
          sr(3)=  -ps(1,3)*(x0(2)-x0(0))*(r0(2)+r0(0))
     &            -ps(2,3)*(x(3)-x0(2)) *(r(3)+r0(2))
     &            -ps(3,3)*(x0(3)-x(3)) *(r0(3)+r(3))
     &            -ps(4,3)*(x0(0)-x0(3))*(r0(0)+r0(3))
          omg(3)= -(x0(2)-x0(0))*(r0(2)+r0(0))
     &            -(x(3)-x0(2)) *(r(3)+r0(2))
     &            -(x0(3)-x(3)) *(r0(3)+r(3))
     &            -(x0(0)-x0(3))*(r0(0)+r0(3))
          sr(4)=  -ps(1,4)*(x0(0)-x0(4))*(r0(0)+r0(4))
     &            -ps(2,4)*(x0(3)-x0(0))*(r0(3)+r0(0))
     &            -ps(3,4)*(x(4)-x0(3)) *(r(4)+r0(3))
     &            -ps(4,4)*(x0(4)-x(4)) *(r0(4)+r(4)) 
          omg(4)= -(x0(0)-x0(4))*(r0(0)+r0(4))
     &            -(x0(3)-x0(0))*(r0(3)+r0(0))
     &            -(x(4)-x0(3)) *(r(4)+r0(3))
     &            -(x0(4)-x(4)) *(r0(4)+r(4)) 
c
          do k=1,4
             if(vol0(k).lt.0d0) then
                vol0(k)=-vol0(k)
                sx(k)=-sx(k)
                sr(k)=-sr(k)
                omg(k)=-omg(k)
             endif
           enddo
c
           if(dden.le.epsilon) cycle
c
           sp(i,j,0)=sp(i,j,0)+sx(1)/(dden*vol0(1))
           sp(i,j,1)=sp(i,j,1)+(sr(1)-omg(1)*pres)/(dden*vol0(1))
           sp(i+1,j,0)=sp(i+1,j,0)+sx(2)/(dden*vol0(2))
           sp(i+1,j,1)=sp(i+1,j,1)+(sr(2)-omg(2)*pres)/(dden*vol0(2))
           sp(i+1,j+1,0)=sp(i+1,j+1,0)+sx(3)/(dden*vol0(3))
           sp(i+1,j+1,1)=sp(i+1,j+1,1)
     1                            +(sr(3)-omg(3)*pres)/(dden*vol0(3))
           sp(i,j+1,0)=sp(i,j+1,0)+sx(4)/(dden*vol0(4))
           sp(i,j+1,1)=sp(i,j+1,1)+(sr(4)-omg(4)*pres)/(dden*vol0(4))
      enddo
      enddo
c
      return
      end 
