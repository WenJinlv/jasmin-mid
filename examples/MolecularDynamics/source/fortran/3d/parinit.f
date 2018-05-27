
c
c************************************************************************
c  初始化数据片.
c************************************************************************
c
      subroutine parinit(xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2, 
     &  gcw0,gcw1, gcw2, 
     &  max_num, real_num,
     &  dbl_attr, dbl_depth, int_attr, int_depth,
     &  index_particles, 
     &  den_gcw0, den_gcw1, den_gcw2,
     &  num_on_cell, is_compute_load)
c
        implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 3, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY,IRZ个分量分别存储x,y,z方向的几何位置.
c 第IVX,IVY,IVZ个分量分别存储x,y,z方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2, IRZ = 3 )
      parameter(IVX = 4, IVY = 5, IVZ = 6 )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY, ICZ)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ICZ,ITYPE, INUM
      parameter(ICX=1,ICY=2,ICZ=3, ITYPE=4, INUM=5)

c
c 计算区域的全局信息. 
c glob_xdx: 单元长度
c glob_xlo: 计算区域左下角坐标
c glob_xhi: 计算区域右上角坐标
c glob_box: 计算区域的全局索引空间
c glob_shift: 周期方向上的单元数.
c period_len: 周期方向上的长度.
c
      double precision glob_xdx(NDIM),glob_xlo(NDIM)
     >  ,glob_xhi(NDIM),glob_xlen(NDIM),period_len(NDIM)
      integer glob_box(2,NDIM), glob_shift(NDIM),myrank
      common/globalData/glob_xdx,glob_xlo,
     >   glob_xhi,glob_xlen,period_len
      common/globalData2/glob_box,glob_shift,myrank
c
c 金属铜的晶格参数
c
      double precision lattice_param(5), lattice_position(NDIM,4), 
     >  lattice_length, cut_length 
      integer num_lattice
      common/CU_param/lattice_param, lattice_position,  
     >  lattice_length, cut_length, num_lattice 
c
      integer ifirst0,ilast0,ifirst1,ilast1, ifirst2,ilast2
      integer gcw0, gcw1, gcw2 
      integer max_num, real_num
      integer den_gcw0, den_gcw1, den_gcw2
      integer dbl_depth, int_depth
      integer is_compute_load
      double precision
     &     xdx(NDIM),xlo(NDIM),xhi(NDIM)
      double precision dbl_attr(0:max_num-1,dbl_depth)

      integer int_attr(0:max_num-1,int_depth)
c
      double precision num_on_cell(
     >  ifirst0-den_gcw0:ilast0+den_gcw0,
     &          ifirst1-den_gcw1:ilast1+den_gcw1,
     &          ifirst2-den_gcw2:ilast2+den_gcw2 )  
c
      integer index_particles(
     >  ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,
     &          ifirst2-gcw2:ilast2+gcw2,INDEX_DEPTH)  
c
c 物理模型参数
c
      double precision triangle_vertex, triangle_angle, 
     >  triangle_depth,tang_value, yycenter,zzcenter
      common/triange/triangle_vertex, triangle_angle, 
     >  triangle_depth,tang_value, yycenter,zzcenter 
      double precision obj_position(2,2),obj_velocity(2)
      integer num_obj, num_atom_in_period(NDIM)
      double precision yyy(2), zzz(2)
      common/Model/obj_position,obj_velocity,yyy,zzz,
     >  num_obj,num_atom_in_period  
c
c 局部变量 
c
      double precision  pr(NDIM,2),pv(NDIM),ra(NDIM), 
     >  omin(NDIM),omax(NDIM), ost(NDIM),err,yy8, 
     >  ytmp,ymin,ymax,step,step9,tmp(NDIM)
      integer i,j,k,io,na,ip,ni,ic,jc,kc,ihead,id,ia
      integer num_atom(NDIM), num_tmp, num_p
      double precision dxdx(NDIM), xyz_start
c
c  begin run
c
        xyz_start = xdx(1)*0.001

        do i=1,NDIM
          dxdx(i) = 1./xdx(i)
        enddo

      ip = 0
      err=EPISLON
C
      do 10 io = 1, num_obj

        pr(1,1) = obj_position(1,io)
        pr(2,1) = obj_position(2,io)
        pv(1) = obj_velocity(io)
        pv(2) = 0.0d0
        pv(3) = 0.0d0
c
        if(pr(2,1).lt.xlo(1).or.pr(1,1).gt.xhi(1)) goto 999
c
        if(pr(1,1).lt.xlo(1))then
          omin(1) = xlo(1)
        else
          omin(1) = pr(1,1)
        endif
c
        if(pr(2,1).gt.xhi(1))then
          omax(1) = xhi(1)  
        else
          omax(1) = pr(2,1) 
        endif
c
        do i=2,NDIM
          omin(i) = xlo(i)
          omax(i) = xhi(i)
        enddo
c
c  compute initiate point in block
c
        do i=1,NDIM
          omin(i) = omin(i) - lattice_length 
          omax(i) = omax(i) + lattice_length 
          k = omin(i)/lattice_length
          ost(i)=k*lattice_length
          num_atom(i) = (omax(i)-ost(i))/lattice_length+2
        enddo
        step = lattice_length
        step9 = 2./step
        do 21 k=0,num_atom(3)
          tmp(3)=ost(3) + step*k + xyz_start
        do 22 j=0,num_atom(2)
          tmp(2)=ost(2) + step*j + xyz_start
        do 23 i=0,num_atom(1)
          tmp(1)=ost(1) + step*i + xyz_start
        do 20 ni=1,1
          do id=1,NDIM
            ra(id) =tmp(id) + lattice_position(id,ni) 
          enddo 
          if( io.eq.1.and.(ra(2).lt.yyy(1).or.ra(2).gt.yyy(2)) )goto 20
c          if( io.eq.1.and.(ra(3).lt.zzz(1).or.ra(3).gt.zzz(2)) )goto 20
 

          if(    ra(1).ge.xlo(1).and.ra(1).lt.xhi(1)
     >      .and.ra(2).ge.xlo(2).and.ra(2).lt.xhi(2)
     >      .and.ra(3).ge.xlo(3).and.ra(3).lt.xhi(3)
     >      .and.tmp(2).lt.period_len(2).and.tmp(3).lt.period_len(3)
     >      .and.ra(1).ge.pr(1,1).and.ra(1).le.pr(2,1))then
c
c  去掉三角形槽中的粒子
c
            if(ra(1).gt.triangle_vertex)then
              yy8 = (ra(1)-triangle_vertex)* tang_value
              ymin = yycenter - yy8
              ymax = yycenter + yy8
              if(ra(2).gt.ymin.and.ra(2).le.ymax)goto 20
            endif                   
C
C record the number of particles in cells.
C
            if(is_compute_load.eq.1)then 
c
              ic = int( ra(1)* dxdx(1) )
              jc = int( ra(2)* dxdx(2) )
              kc = int( ra(3)* dxdx(3) )
             if(ic.ge.ifirst0.and.ic.le.ilast0.and. 
     >         jc.ge.ifirst1.and.jc.le.ilast1.and.
     >         kc.ge.ifirst2.and.kc.le.ilast2)then
c
               num_on_cell(ic,jc,kc) = num_on_cell(ic,jc,kc) + 1.0             
              endif
            else
              dbl_attr(ip, IRX) = ra(1)
              dbl_attr(ip, IRY) = ra(2)
              dbl_attr(ip, IRZ) = ra(3)
              dbl_attr(ip, IVX) = pv(1)
              dbl_attr(ip, IVY) = pv(2)
              dbl_attr(ip, IVZ) = pv(3)
              int_attr(ip, ITYPE) = 1
              int_attr(ip, INUM) = ip
            endif
            ip=ip+1
c
          endif
20      continue
23      continue
22      continue
21      continue
999   continue
c  循环设置两个物体
10    continue

c
      real_num = ip
      if(is_compute_load.eq.1)return 
c
c  计算粒子所属单元的索引号.
c
      call computeCellOfParticles(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,glob_box,glob_shift,
     >  xdx, real_num,
     >  dbl_attr(0,IRX),dbl_attr(0,IRY),dbl_attr(0,IRZ),
     >  int_attr(0,ICX),int_attr(0,ICY),int_attr(0,ICZ) ) 
c

      return
      end   
c
c***********************************************************************
c  设置全局信息
c***********************************************************************     
      subroutine setglobalparam(
     &  myrank8,shift8,xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  left_obj_vel,left_obj_pos, 
     &  right_obj_vel,right_obj_pos,
     &  vertex, angle) 
c
        implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 3, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY,IRZ个分量分别存储x,y,z方向的几何位置.
c 第IVX,IVY,IVZ个分量分别存储x,y,z方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2, IRZ = 3 )
      parameter(IVX = 4, IVY = 5, IVZ = 6 )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY, ICZ)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ICZ,ITYPE, INUM
      parameter(ICX=1,ICY=2,ICZ=3, ITYPE=4, INUM=5)

c
c 计算区域的全局信息. 
c glob_xdx: 单元长度
c glob_xlo: 计算区域左下角坐标
c glob_xhi: 计算区域右上角坐标
c glob_box: 计算区域的全局索引空间
c glob_shift: 周期方向上的单元数.
c period_len: 周期方向上的长度.
c
      double precision glob_xdx(NDIM),glob_xlo(NDIM)
     >  ,glob_xhi(NDIM),glob_xlen(NDIM),period_len(NDIM)
      integer glob_box(2,NDIM), glob_shift(NDIM),myrank
      common/globalData/glob_xdx,glob_xlo,
     >   glob_xhi,glob_xlen,period_len
      common/globalData2/glob_box,glob_shift,myrank
c
c 金属铜的晶格参数
c
      double precision lattice_param(5), lattice_position(NDIM,4), 
     >  lattice_length, cut_length 
      integer num_lattice
      common/CU_param/lattice_param, lattice_position,  
     >  lattice_length, cut_length, num_lattice 
c 
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
c
      double precision
     &     xdx(NDIM),xlo(NDIM),xhi(NDIM)
      integer shift8(NDIM),myrank8
      double precision vertex, angle, depth
      double precision left_obj_vel,left_obj_pos(2) 
      double precision right_obj_vel,right_obj_pos(2) 
c
c 物理模型参数
c
      double precision triangle_vertex, triangle_angle, 
     >    triangle_depth,tang_value, yycenter,zzcenter
      common/triange/triangle_vertex, triangle_angle, 
     >    triangle_depth,tang_value, yycenter,zzcenter 
      double precision obj_position(2,2),obj_velocity(2)
      integer num_obj, num_atom_in_period(NDIM)
      double precision yyy(2), zzz(2)
      common/Model/obj_position,obj_velocity, yyy, zzz,
     >    num_obj,num_atom_in_period  
c 
c 局部变量
c
      integer i,j,k,ii,ic,jc
c
c begin   
c
      myrank = myrank8
      do i=1, NDIM
        glob_xdx(i) = xdx(i)         
        glob_xlo(i) = xlo(i)         
        glob_xhi(i) = xhi(i)         
        glob_shift(i) = shift8(i)         
      enddo
c global domain bounding box
      glob_box(1,1) = ifirst0
      glob_box(2,1) = ilast0
      glob_box(1,2) = ifirst1
      glob_box(2,2) = ilast1
      glob_box(1,3) = ifirst2
      glob_box(2,3) = ilast2

      do i=1, NDIM
        glob_xlen(i) = xhi(i) - xlo(i)
      enddo

      depth = 0.5*(left_obj_pos(2)-left_obj_pos(1))
      yyy(1)= glob_xlen(2) * (0.5 - depth)
      yyy(2)= glob_xlen(2) * (0.5 + depth)
      zzz(1)= glob_xlen(3) * (0.5 - depth)
      zzz(2)= glob_xlen(3) * (0.5 + depth)
      if( (zzz(2) - zzz(1) )< 2*xdx(3) )then
         zzz(1) = 0
         zzz(2) = glob_xlen(3)
      endif
c
c 设置物理模型
c
      num_obj = 2
c  左边的物体
      obj_velocity(1) = left_obj_vel 
      obj_position(1,1) = left_obj_pos(1) * glob_xlen(1)
      obj_position(2,1) = left_obj_pos(2) * glob_xlen(1)
c  右边的物体
      obj_velocity(2) = right_obj_vel 
      obj_position(1,2) = right_obj_pos(1) * glob_xlen(1)
      obj_position(2,2) = right_obj_pos(2) * glob_xlen(1)
c  三角形槽
      triangle_vertex = vertex 
      triangle_angle  = angle
      triangle_depth  = 0.1
      triangle_vertex = triangle_vertex * glob_xlen(1)
      triangle_depth = triangle_depth * glob_xlen(1)
      tang_value = tan(0.5*triangle_angle/180*PI)
      yycenter = glob_xlen(2)*0.5
      zzcenter = glob_xlen(3)*0.5
c
c   金属铜的晶格参数. 
c
      lattice_param(1) = 1.0550d0
      lattice_param(2) = 54.9300d0
      lattice_param(3) = 1.3588d0
      lattice_param(4) = 2.8660d0
      lattice_length   = 3.615D0
      cut_length       = lattice_length*3.d0 + xdx(1)*0.01
c
c  三维面心立方晶胞中的原子位置.
c
       num_lattice = 4
       do i=1,num_lattice
         do j=1,NDIM
           lattice_position(j,i)=0.0d0
         enddo
       enddo
       lattice_position(2,2)=0.5d0*lattice_length
       lattice_position(3,2)=lattice_position(2,2)
       lattice_position(1,3)=lattice_position(2,2)
       lattice_position(3,3)=lattice_position(2,2)
       lattice_position(1,4)=lattice_position(2,2)
       lattice_position(2,4)=lattice_position(2,2)
c
c 周期方向上, 原子个数和周期长度.
c
      do i=1,NDIM
        ii = glob_xlen(i)/lattice_length
c        period_len(i) = lattice_length * ii
        period_len(i) = glob_xlen(i)
        num_atom_in_period(i) = ii
      enddo
c
c 修正物理模型
c
      do i =1, num_obj
        ii = obj_position(1,i)/lattice_length
        obj_position(1,i) = ii * lattice_length 
        ii = obj_position(2,i)/lattice_length
        obj_position(2,i) = ii * lattice_length 
      enddo
c
      if(myrank.eq.0.and.DEBUG_CHECK)then
        open(88,file='record.dat')
        write(88,*)'---------------------------------------------'
        write(88,*)'----  computational domain parameter:   -----'
        write(88,*)'physical domain: ', (glob_xlen(i), i=1,NDIM)
        write(88,*)'logical domain:  ', 
     >    (glob_box(1,i), glob_box(2,i),i=1,NDIM)
        write(88,*)'cell length:     ', (glob_xdx(i),i=1,NDIM)
        write(88,*)'periodic shift:  ', (glob_shift(i),i=1,NDIM)
        write(88,*)'---------------------------------------------'
        write(88,*)'----  physical model parameter:   -----'
        write(88,*)'lattice length: ', lattice_length
        write(88,*)'cut off length: ', cut_length
c
c  检测：单元边长大于截断距离.
c
        do i=1,NDIM
          if(myrank.eq.0.and.xdx(i).le.cut_length)then
            write(*,*)'Error!!! cell length <= cut off length' 
            write(*,*)'cell length >', cut_length+epislon 
            write(*,*)(xdx(k),k=1,NDIM), cut_length
            stop 777
          endif
        enddo
        write(88,*)'periodic length: ', (period_len(i),i=1,NDIM)
        write(88,*)'Max atom in XYZ: ', 
     >   (num_atom_in_period(i),i=1,NDIM)
        write(88,*)'---------------------------------------------'
        write(88,*)'object position and velocity: ', num_obj 
        do i =1, num_obj
          write(88,*)obj_position(1,i),
     >      obj_position(2,i),obj_velocity(i)
        enddo
        close(88)
c
      endif
c
      return
      end
      
c*********************************************************************
c  计算网格片单元上的负载权重
c*********************************************************************
      subroutine computeload(ttt, xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2, 
     &  gcw0,gcw1,gcw2, lgcw0,lgcw1,lgcw2, 
     &  index_particles, load_weight)
c
        implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 3, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY,IRZ个分量分别存储x,y,z方向的几何位置.
c 第IVX,IVY,IVZ个分量分别存储x,y,z方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2, IRZ = 3 )
      parameter(IVX = 4, IVY = 5, IVZ = 6 )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY, ICZ)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ICZ,ITYPE, INUM
      parameter(ICX=1,ICY=2,ICZ=3, ITYPE=4, INUM=5)

c
c 计算区域的全局信息. 
c glob_xdx: 单元长度
c glob_xlo: 计算区域左下角坐标
c glob_xhi: 计算区域右上角坐标
c glob_box: 计算区域的全局索引空间
c glob_shift: 周期方向上的单元数.
c period_len: 周期方向上的长度.
c
      double precision glob_xdx(NDIM),glob_xlo(NDIM)
     >  ,glob_xhi(NDIM),glob_xlen(NDIM),period_len(NDIM)
      integer glob_box(2,NDIM), glob_shift(NDIM),myrank
      common/globalData/glob_xdx,glob_xlo,
     >   glob_xhi,glob_xlen,period_len
      common/globalData2/glob_box,glob_shift,myrank
c
c 金属铜的晶格参数
c
      double precision lattice_param(5), lattice_position(NDIM,4), 
     >  lattice_length, cut_length 
      integer num_lattice
      common/CU_param/lattice_param, lattice_position,  
     >  lattice_length, cut_length, num_lattice 
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0, gcw1,gcw2,lgcw0,lgcw1,lgcw2
      double precision ttt, xdx(NDIM),xlo(NDIM),xhi(NDIM)

      double precision  
     >  load_weight(ifirst0-lgcw0:ilast0+lgcw0,
     &          ifirst1-lgcw1:ilast1+lgcw1,
     &          ifirst2-lgcw2:ilast2+lgcw2,1)  
C
      integer  
     >  index_particles(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,
     &          ifirst2-gcw2:ilast2+gcw2,INDEX_DEPTH)  

c
c 局部变量 
c
      integer i0, i1,i2

      do i2= ifirst2,ilast2
      do i1= ifirst1,ilast1
      do i0= ifirst0,ilast0
        load_weight(i0,i1,i2,1) = index_particles(i0,i1,i2,1)*1.0
      enddo
      enddo
      enddo
c
      return
      end   
