define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c************************************************************************
c  计算粒子受力，更新粒子位置.
c************************************************************************
      subroutine updateparticles(dt8,xdx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  max_num, inner_num,
     &  dbl_attr, dbl_depth, 
     &  int_attr, int_depth, 
     &  index_particles_by_cell,
     &  pf)
c
      implicit none
include(FORTDIR/../3d/model_parameter.i)dnl
c
      integer ifirst0,ilast0,ifirst1,ilast1, ifirst2,ilast2
      integer gcw0, gcw1, gcw2 
      integer max_num, inner_num
      integer dbl_depth, int_depth 
      logical non_ghost_cell
      REAL
     &     dt8, xdx(NDIM), xlo(NDIM), xhi(NDIM)

      REAL 
     > dbl_attr(0:max_num-1,dbl_depth)
     > , pf(0:inner_num-1,NDIM)

       integer index_particles_by_cell
     >   (CELL3dVECG(ifirst,ilast,gcw), INDEX_DEPTH)
c
      integer int_attr(0:max_num-1,int_depth)
c
c 局部变量
c
      double precision rrcut2,rr_ref,rr,fix,f_ref,ftol
      double precision fx1,fx2,fx3,drx,dry,drz,f 
      double precision da1,da2,dmass 
      integer mtol,ia,ic,jc,kc,ipos,ic2,jc2,kc2
      integer ipa2,ia2,ipos2,ipa,num_tmp
      integer na,na2,na_end,na_begin,na2_end,na2_begin
c
c 计算参考距离处的受力. 在参考距离外直至截断距离, 受力线性下降为零.
c
      rrcut2 = cut_length * cut_length 
      rr_ref = rrcut2 * 0.9 
      rr     = dsqrt(rr_ref)
      fix    = dexp(-lattice_param(3)*(rr-lattice_param(4))) 
      f_ref  = -2*lattice_param(2)*lattice_param(3)*fix*(fix-1.0d0)
c
c  计算粒子受力
c
      do ia = 0,inner_num-1
        pf(ia, 1) = 0.0d0
        pf(ia, 2) = 0.0d0
        pf(ia, 3) = 0.0d0
      enddo
c
c 遍历网格片中的计算单元: 内部单元和影像单元, 不包括影像单元. 
c
      do 100 kc = ifirst2, ilast2
      do 100 jc = ifirst1, ilast1
      do 100 ic = ifirst0, ilast0
        na   = index_particles_by_cell(ic,jc,kc,1)
        ipos = index_particles_by_cell(ic,jc,kc,2)-1
        if(na.lt.1)goto 100
c
c 循环: (ic,jc,kc)的相邻单元.  
c
      do 200 kc2 = kc-1, kc+1
      do 200 jc2 = jc-1, jc+1
      do 200 ic2 = ic-1, ic+1
c
c  基于作用力与反作用力计算粒子受力
c  如果相邻单元均是内部单元, 则搜索周围13个单元和本单元内的粒子: 
c  (ic-1,jc+1,kc+1), (ic,jc+1,kc+1), (ic+1,jc+1,kc+1)
c  (ic-1,jc,  kc+1),   (ic,jc,kc+1),   (ic+1,jc,kc+1)
c  (ic-1,jc-1,kc+1), (ic,jc-1,kc+1), (ic+1,jc-1,kc+1)
c  (ic-1,jc+1,kc),   (ic,jc+1,kc),   (ic+1,jc+1,kc)
c                    (ic,jc,kc)  ,   (ic+1,jc,kc), 
c  如果相邻单元是影像单元, 则必须计算, 但不将力累加到影像单元中的粒子. 
c
        non_ghost_cell=.true.
        if(kc2.lt.ifirst2.or.kc2.gt.ilast2.or.
     >    jc2.lt.ifirst1.or.jc2.gt.ilast1.or.
     >    ic2.lt.ifirst0.or.ic2.gt.ilast0)
     >    non_ghost_cell=.false.
          if(non_ghost_cell)then
            if(kc2.eq.kc-1)then
            goto 200
          else if(kc2.eq.kc)then 
            if(jc2.eq.jc-1.or.(jc2.eq.jc.and.ic2.eq.ic-1))goto 200
          endif
        endif
c
        na2 = index_particles_by_cell(ic2,jc2,kc2,1)
        ipos2 = index_particles_by_cell(ic2,jc2,kc2,2)-1
        if(ic2.eq.ic.and.jc2.eq.jc.and.kc2.eq.kc)then   
          na_begin = 1
          na_end = na - 1
        else
          na_begin = 1
          na_end = na
        endif
c
        do 300 ia = na_begin, na_end 
          if(ic2.eq.ic.and.jc2.eq.jc.and.kc2.eq.kc)then   
            na2_begin = ia+1
          else
            na2_begin = 1
          endif
          na2_end = na2
c
          ipa = ipos + ia
          fx1 = 0.0d0
          fx2 = 0.0d0
          fx3 = 0.0d0
          do 400 ia2 = na2_begin, na2_end 
            ipa2 = ipos2 + ia2
            drx = dbl_attr(ipa2,IRX)-dbl_attr(ipa,IRX)
            dry = dbl_attr(ipa2,IRY)-dbl_attr(ipa,IRY)
            drz = dbl_attr(ipa2,IRZ)-dbl_attr(ipa,IRZ)
            rr  = drx*drx+dry*dry+drz*drz
            if(rr.gt.rrcut2) goto 400 
c
c 调试代码: 如果两个粒子的位置重叠, 则终止程序.
c
            if(DEBUG_CHECK)then 
            if(rr.lt.EPISLON)then
              write(88,*)' Error message!!!'
              write(88,*)ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
              write(88,*)ipa,ipa2,ipos,ipos2
              write(88,*)dbl_attr(ipa,IRX),dbl_attr(ipa,IRY)
     >                  ,dbl_attr(ipa,IRZ),ic,jc
              write(88,*)dbl_attr(ipa2,IRX),dbl_attr(ipa2,IRY)
     >                  ,dbl_attr(ipa2,IRZ),ic2,jc2
              close(88)
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
              f = f_ref * (rrcut2-rr)/(rrcut2-rr_ref)
            endif

            f   = f/rr
            fx1 = fx1 + f*drx
            fx2 = fx2 + f*dry
            fx3 = fx3 + f*drz

            if(non_ghost_cell)then
              pf(ipa2,1) = pf(ipa2,1) - f*drx
              pf(ipa2,2) = pf(ipa2,2) - f*dry
              pf(ipa2,3) = pf(ipa2,3) - f*drz
            endif
400       continue
          pf(ipa,1) = pf(ipa,1) + fx1
          pf(ipa,2) = pf(ipa,2) + fx2
          pf(ipa,3) = pf(ipa,3) + fx3
300     continue
200   continue
100   continue
c
c 计算粒子速度, 并更新粒子位置.
c
      dmass = 1.0d0/lattice_param(1) * dt8
      do 500 ia = 0, inner_num-1
        dbl_attr(ia,IVX) = dbl_attr(ia,IVX) + pf(ia,1)*dmass
        dbl_attr(ia,IVY) = dbl_attr(ia,IVY) + pf(ia,2)*dmass
        dbl_attr(ia,IVZ) = dbl_attr(ia,IVZ) + pf(ia,3)*dmass
500   continue
      do 600 ia = 0, inner_num-1
        dbl_attr(ia,IRX) = dbl_attr(ia,IRX) + dbl_attr(ia,IVX) *dt8
        dbl_attr(ia,IRY) = dbl_attr(ia,IRY) + dbl_attr(ia,IVY) *dt8
        dbl_attr(ia,IRZ) = dbl_attr(ia,IRZ) + dbl_attr(ia,IVZ) *dt8
600   continue
c
      call computeCellOfParticles(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,glob_box,glob_shift,
     >  xdx, inner_num,
     >  dbl_attr(0,IRX),dbl_attr(0,IRY),dbl_attr(0,IRZ),
     >  int_attr(0,ICX),int_attr(0,ICY),int_attr(0,ICZ) ) 
c
      return
      end   
c
c***********************************************************************
c  计算粒子所属单元的索引号 
c***********************************************************************     
c
      subroutine computeCellOfParticles(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2, glob_box,glob_shift,
     >  xdx, real_num,
     >  pos_x, pos_y, pos_z,cell_x, cell_y, cell_z)
c         
       integer  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
       integer  gcw0,gcw1,gcw2,glob_box(2,3),glob_shift(3)
       integer real_num
       double precision xdx(3), ddx, ddy, ddz
       double precision pos_x(real_num)          
       double precision pos_y(real_num)          
       double precision pos_z(real_num)          
       integer cell_x(real_num)          
       integer cell_y(real_num)          
       integer cell_z(real_num)          
C
       integer ia
C         
       ddx = 1./xdx(1)     
       ddy = 1./xdx(2)     
       ddz = 1./xdx(3)     
       do ia = 1, real_num
          cell_x(ia) = int( pos_x(ia) * ddx +100) -100
          cell_y(ia) = int( pos_y(ia) * ddy +100) -100
          cell_z(ia) = int( pos_z(ia) * ddz +100) -100
       enddo                 
C
C tag particles out of computational domain
C
C x direction
C
      if(glob_shift(1).eq.0.and.ifirst0.eq.glob_box(1,1))then
        do ia = 1, real_num
          if( cell_x(ia) .lt. glob_box(1,1) )
     >       cell_x(ia) = cell_x(ia) -gcw0 
        enddo                 
      endif                 
      if(glob_shift(1).eq.0.and.ilast0.eq.glob_box(2,1))then
        do ia = 1, real_num
          if( cell_x(ia) .gt. glob_box(2,1) )
     >       cell_x(ia) = cell_x(ia) + gcw0 
        enddo                 
      endif                 
C
C y direction
C
      if(glob_shift(2).eq.0.and.ifirst1.eq.glob_box(1,2))then
        do ia = 1, real_num
          if( cell_y(ia) .lt. glob_box(1,2) )
     >       cell_y(ia) = cell_y(ia) -gcw1 
        enddo                 
      endif                 
      if(glob_shift(2).eq.0.and.ilast1.eq.glob_box(2,2))then
        do ia = 1, real_num
          if( cell_y(ia) .gt. glob_box(2,2) )
     >       cell_y(ia) = cell_y(ia) + gcw1 
        enddo                 
      endif                 
C
C z direction
C
      if(glob_shift(3).eq.0.and.ifirst2.eq.glob_box(1,3))then
        do ia = 1, real_num
          if( cell_z(ia) .lt. glob_box(1,3) )
     >       cell_z(ia) = cell_z(ia) -gcw2 
        enddo                 
      endif                 
      if(glob_shift(3).eq.0.and.ilast2.eq.glob_box(2,3))then
        do ia = 1, real_num
          if( cell_z(ia) .gt. glob_box(2,3) )
     >       cell_z(ia) = cell_z(ia) + gcw2 
        enddo                 
      endif                 

       return
       end

c
c************************************************************************
c  计算网格片上的粒子总动量(x方向).
c************************************************************************
C
      subroutine computemomentum(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2, 
     &  gcw0,gcw1,gcw2, psize,dbl_depth, 
     &  dbl_attr, index_particles_by_cell,
     &  momentum)
c
      implicit none
include(FORTDIR/../3d/model_parameter.i)dnl
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0, gcw1, gcw2, psize, dbl_depth
      double precision dbl_attr(dbl_depth,0:psize) 
       integer index_particles_by_cell
     >   (CELL3dVECG(ifirst,ilast,gcw), INDEX_DEPTH)
      double precision momentum
c
      integer ic,jc,kc,na,ibase,ia,ip

c
c 计算x方向的动量
c
      momentum = 0.0d0
      do kc = ifirst2, ilast2
      do jc = ifirst1, ilast1
      do ic = ifirst0, ilast0
        na = index_particles_by_cell(ic,jc,kc,1)
        ibase = index_particles_by_cell(ic,jc,kc,2) -1
        do ia = 1, na 
          ip = ibase + ia 
          momentum = momentum + dbl_attr(IVX,ip)
        enddo
      enddo
      enddo
      enddo
c
      return
      end
c
      subroutine computedensity(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2, 
     &  pgcw0, pgcw1, pgcw2, index_particles,
     &  dgcw0, dgcw1, dgcw2, particle_density)
c
        implicit none
include(FORTDIR/../3d/model_parameter.i)dnl
c
      integer ifirst0,ilast0,ifirst1,ilast1, ifirst2,ilast2
      integer dgcw0, dgcw1, dgcw2 
      integer pgcw0, pgcw1, pgcw2 
       REAL   
     >  particle_density(CELL3dVECG(ifirst,ilast,dgcw) )
      integer  index_particles(
     >  CELL3dVECG(ifirst,ilast,pgcw),INDEX_DEPTH)  
c
      integer kc, jc, ic, num_tmp 
c
c 计算单元上的粒子密度
c
      do kc = ifirst2, ilast2
      do jc = ifirst1, ilast1
      do ic = ifirst0, ilast0
        num_tmp = index_particles(ic,jc,kc,1)
        particle_density(ic,jc,kc) = num_tmp*1.0 
      enddo
      enddo
      enddo
C
      return
      end

