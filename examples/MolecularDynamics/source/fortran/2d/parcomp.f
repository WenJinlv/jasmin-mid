









c
c************************************************************************
c  计算粒子受力，更新粒子位置. 
c************************************************************************
      subroutine updateparticles(dt8,xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1, 
     &  gcw0,gcw1, 
     &  max_num, inner_num,
     &  dbl_attr, dbl_depth, 
     &  int_attr, int_depth, 
     &  index_particles_by_cell,
     &  pf)
c
      implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 2, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY个分量分别存储x,y方向的几何位置.
c 第IVX,IVY个分量分别存储x,y方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2  )
      parameter(IVX = 3, IVY = 4  )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ITYPE, INUM
      parameter(ICX=1,ICY=2, ITYPE=3, INUM=4)

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
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0, gcw1
      integer max_num, inner_num
      integer dbl_depth, int_depth 
      logical non_ghost_cell
      double precision
     &     dt8, xdx(NDIM), xlo(NDIM), xhi(NDIM)

      double precision 
     > dbl_attr(0:max_num-1,dbl_depth)
     > , pf(0:inner_num-1,NDIM)

       integer index_particles_by_cell
     >   (ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1, INDEX_DEPTH)
c
      integer int_attr(0:max_num-1,int_depth)
c
c 局部变量
c
      double precision rrcut2,rr_ref,rr,fix,f_ref,ftol
      double precision fx1,fx2,drx,dry,drz,f 
      double precision da1,da2,dmass 
      integer mtol,ia,ic,jc,ipos,ic2,jc2
      integer ipa2,ia2,ipos2,ipa,num_tmp
      integer na,na2,na_end,na_begin,na2_end,na2_begin
c
c 计算参考距离处的受力. 在参考距离直至截断距离处, 受力线性下降为零.
c
      rrcut2 = cut_length * cut_length 
      rr_ref = rrcut2 * 0.9 
      rr     = dsqrt(rr_ref)
      f      = (lattice_param(4)-rr)*lattice_param(3)
      fix    = dexp(f) 

      f_ref  = -2*lattice_param(2)*lattice_param(3)*fix*(fix-1.0d0)
c
c  计算粒子受力
c
      ftol = 0.0d0
      mtol = 0
      do ia = 0,inner_num-1
        pf(ia,1) = 0.0d0
        pf(ia,2) = 0.0d0
      enddo
c
c 遍历网格片中的内部单元. 
c
      do 100 jc = ifirst1, ilast1
      do 100 ic = ifirst0, ilast0
        na   = index_particles_by_cell(ic,jc,1)
        ipos = index_particles_by_cell(ic,jc,2)-1
c
c 循环: (ic,jc)的相邻单元.  
c
      do 200 jc2 = jc-1, jc+1
      do 200 ic2 = ic-1, ic+1
c
c  基于作用力与反作用力定律计算粒子受力
c  如果相邻单元均是内部单元, 则计算(ic-1,jc+1), (ic,jc+1), 
c  (ic+1,jc+1), (ic+1,jc), (ic,jc)五个单元.
c  如果相邻单元是影像单元, 则必须计算, 但不将力累加到影像单元中的粒子. 
c
        non_ghost_cell = .true.
        if(jc2.eq.ifirst1-1.or.jc2.eq.ilast1+1.or.
     >    ic2.eq.ifirst0-1.or.ic2.eq.ilast0+1)non_ghost_cell=.false.
        if(non_ghost_cell)then
          if(jc2.eq.jc-1.or.(jc2.eq.jc.and.ic2.eq.ic-1))goto 200
        endif
c
        na2   = index_particles_by_cell(ic2,jc2,1)
        ipos2 = index_particles_by_cell(ic2,jc2,2)-1
        if(ic2.eq.ic.and.jc2.eq.jc)then   
c 本单元内部粒子相互作用
          na_begin = 1
          na_end   = na - 1
        else
c 两个单元之间的粒子相互作用
          na_begin = 1
          na_end   = na
        endif
c 
        do 300 ia = na_begin, na_end 
          if(ic2.eq.ic.and.jc2.eq.jc)then   
c 本单元内部粒子相互作用
            na2_begin = ia+1
          else
c 两个单元之间的粒子相互作用
            na2_begin = 1
          endif
          na2_end = na2
c
          ipa = ipos + ia
          fx1 = 0.0d0
          fx2 = 0.0d0
          do 400 ia2 = na2_begin, na2_end 
            ipa2 = ipos2 + ia2
            drx = dbl_attr(ipa2,IRX)-dbl_attr(ipa,IRX)
            dry = dbl_attr(ipa2,IRY)-dbl_attr(ipa,IRY)
            rr  = drx*drx+dry*dry
            if(rr.gt.rrcut2) goto 400 
c
c 调试代码: 如果两个粒子的位置重叠, 则终止程序.
c
            if(DEBUG_CHECK)then 
            if(rr.lt.EPISLON)then
              write(*,*)' Error message in parcomp.f'
              write(*,*)' the distance of two atom is too close'
              write(*,*)ifirst0,ilast0,ifirst1,ilast1
              write(*,*)dbl_attr(ipa,irx),dbl_attr(ipa,iry),ic,jc
              write(*,*)dbl_attr(ipa2,irx),dbl_attr(ipa2,iry),ic2,jc2
              stop 777
            endif
            endif
c
            if(rr.lt.rr_ref)then
c
c 根据Morse势计算粒子受力.
c
              rr  = dsqrt(rr)
              fix = dexp(-lattice_param(3)*(rr-lattice_param(4))) 
              f   = -2*lattice_param(2)*lattice_param(3)*fix*(fix-1.0d0)
            else
c
c 修正参考距离外的力, 以确保在截断距离处的力为零.
c
              f   = f_ref * (rrcut2-rr)/(rrcut2-rr_ref)
            endif

            f   = f/rr
            fx1 = fx1+f*drx
            fx2 = fx2+f*dry
            if(non_ghost_cell)then
              pf(ipa2,1) = pf(ipa2,1)-f*drx
              pf(ipa2,2) = pf(ipa2,2)-f*dry
            endif
400       continue
          pf(ipa,1) = pf(ipa,1)+fx1
          pf(ipa,2) = pf(ipa,2)+fx2
300     continue
200   continue
100   continue
c
c 计算粒子速度, 并更新粒子位置.
c
      dmass = 1.0d0/lattice_param(1) * dt8
      do 500 ia = 0, inner_num-1
        da1 = pf(ia,1)*dmass
        da2 = pf(ia,2)*dmass
        dbl_attr(ia,IVX) = dbl_attr(ia,IVX) + pf(ia,1)*dmass
        dbl_attr(ia,IVY) = dbl_attr(ia,IVY) + pf(ia,2)*dmass
500   continue
      do 600 ia = 0, inner_num-1
        dbl_attr(ia,IRX) = dbl_attr(ia,IRX) + dbl_attr(ia,IVX) *dt8
        dbl_attr(ia,IRY) = dbl_attr(ia,IRY) + dbl_attr(ia,IVY) *dt8
600   continue
c
      call computeCellOfParticles(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,glob_box,glob_shift,
     >  xdx,inner_num,
     >  dbl_attr(0,IRX),dbl_attr(0,IRY),
     >  int_attr(0,ICX),int_attr(0,ICY) ) 
c
      return
      end   
c***********************************************************************
c  计算粒子所属的单元号
c***********************************************************************     

      subroutine computeCellOfParticles(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,glob_box,glob_shift,
     >  xdx, real_num_of_particles,
     >  pos_x, pos_y, cell_x, cell_y)
C         
       integer  ifirst0,ilast0,ifirst1,ilast1
       integer  gcw0,gcw1,glob_box(2,2),glob_shift(2)
       integer real_num_of_particles
       double precision xdx(2), ddx, ddy
       double precision pos_x(real_num_of_particles)          
       double precision pos_y(real_num_of_particles)          
       integer cell_x(real_num_of_particles)          
       integer cell_y(real_num_of_particles)          
C
       integer ia
C         
       ddx = 1./xdx(1)     
       ddy = 1./xdx(2)     
       do ia = 1, real_num_of_particles
          cell_x(ia) = int( pos_x(ia) * ddx +2) -2
          cell_y(ia) = int( pos_y(ia) * ddy +2) -2 
       enddo                 
C
C tag particles out of computational domain
C
C x direction
C
      if(glob_shift(1).eq.0.and.ifirst0.eq.glob_box(1,1))then
        do ia = 1, real_num_of_particles
          if( cell_x(ia) .lt. glob_box(1,1) )
     >       cell_x(ia) = cell_x(ia) -gcw0 
        enddo                 
      endif                 
      if(glob_shift(1).eq.0.and.ilast0.eq.glob_box(2,1))then
        do ia = 1, real_num_of_particles
          if( cell_x(ia) .gt. glob_box(2,1) )
     >       cell_x(ia) = cell_x(ia) + gcw0 
        enddo                 
      endif                 
C
C y direction
C
      if(glob_shift(2).eq.0.and.ifirst1.eq.glob_box(1,2))then
        do ia = 1, real_num_of_particles
          if( cell_y(ia) .lt. glob_box(1,2) )
     >       cell_y(ia) = cell_y(ia) -gcw1 
        enddo                 
      endif                 
      if(glob_shift(2).eq.0.and.ilast1.eq.glob_box(2,2))then
        do ia = 1, real_num_of_particles
          if( cell_y(ia) .gt. glob_box(2,2) )
     >       cell_y(ia) = cell_y(ia) + gcw1 
        enddo                 
      endif                 
C
       return
       end
c
      subroutine computeload(ttt, xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1, 
     &  gcw0,gcw1, lgcw0,lgcw1, 
     &  index_particles, load_weight)
c
        implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 2, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY个分量分别存储x,y方向的几何位置.
c 第IVX,IVY个分量分别存储x,y方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2  )
      parameter(IVX = 3, IVY = 4  )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ITYPE, INUM
      parameter(ICX=1,ICY=2, ITYPE=3, INUM=4)

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
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0, gcw1,lgcw0,lgcw1
      double precision ttt, xdx(NDIM),xlo(NDIM),xhi(NDIM)

      double precision  
     >  load_weight(ifirst0-lgcw0:ilast0+lgcw0,
     &          ifirst1-lgcw1:ilast1+lgcw1,1)  
C
      integer  
     >  index_particles(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,INDEX_DEPTH)  

c
c 局部变量 
c
      integer i0, i1,i2

      do i1= ifirst1,ilast1
      do i0= ifirst0,ilast0
        load_weight(i0,i1,1) = index_particles(i0,i1,1)*1.0+1.0
      enddo
      enddo
c
      return
      end   
      
      subroutine computedensity(
     &  ifirst0,ilast0,ifirst1,ilast1, 
     &  pgcw0,pgcw1,index_particles,  
     &  dgcw0,dgcw1,particle_density)
c
        implicit none
c
c  常量
c
      logical DEBUG_CHECK
      parameter(DEBUG_CHECK = .false.)

      integer NDIM, INDEX_DEPTH
      parameter(NDIM = 2, INDEX_DEPTH = 2)

      double precision PI,EPISLON
      parameter(PI = 3.1415926)
      parameter(EPISLON = 1.0D-10)
c 
c 粒子的属性按数据类型分为浮点型和整型, 分别存储在 
c 浮点型属性数组dbl_attr和整型属性数组int_attr中.
c
c 浮点型属性数组：存储粒子的6个属性.
c 第IRX,IRY个分量分别存储x,y方向的几何位置.
c 第IVX,IVY个分量分别存储x,y方向的速度.
c
      integer IRX,IRY,IRZ,IVX,IVY,IVZ
      parameter(IRX = 1, IRY = 2  )
      parameter(IVX = 3, IVY = 4  )
c 
c 整型属性数组：存储粒子的整型属性.
c 前NDIM个分量(ICX, ICY)存储粒子所属单元的索引号.
c 第ITYPE个分量存储粒子的物质类型.
c 第INUM个分量存储粒子的编号.
c
      integer ICX,ICY,ITYPE, INUM
      parameter(ICX=1,ICY=2, ITYPE=3, INUM=4)

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
      integer ifirst0,ilast0,ifirst1,ilast1
      integer dgcw0, dgcw1,pgcw0,pgcw1

      double precision  
     >  particle_density(ifirst0-dgcw0:ilast0+dgcw0,
     &          ifirst1-dgcw1:ilast1+dgcw1,1)  
C
      integer  
     >  index_particles(ifirst0-pgcw0:ilast0+pgcw0,
     &          ifirst1-pgcw1:ilast1+pgcw1
     >   ,INDEX_DEPTH)  

c
c 局部变量 
c
      integer i0, i1,i2

      do i1= ifirst1,ilast1
      do i0= ifirst0,ilast0
        particle_density(i0,i1,1) = index_particles(i0,i1,1)*1.0
      enddo
      enddo
c
      return
      end   

