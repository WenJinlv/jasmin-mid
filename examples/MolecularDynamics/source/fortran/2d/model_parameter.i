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
