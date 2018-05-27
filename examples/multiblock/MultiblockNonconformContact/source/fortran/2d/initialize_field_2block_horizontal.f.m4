include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c     
c    本初始化子程序用于包含两个block的sod模型. 两个block垂直相邻分布,
c    两block的交接面垂直于Y轴, 用滑移线刻画.
c     
c***********************************************************************
      subroutine initializefield2blockh(
     &     ifirst0, ilast0,
     &     ifirst1, ilast1,
     &     gcw0, gcw1,
     &     block_id,
     &     d_velocity_in,
     &     velocity,
     &     node_coord)
c***********************************************************************
      implicit none
c***********************************************************************
c input parameter:
c
      integer ifirst0, ilast0, ifirst1, ilast1,
     &        gcw0, gcw1, block_id
      REAL d_velocity_in(1:NDIM),
     &     node_coord(NODE2dVECG(ifirst,ilast,gcw),1:NDIM)
c     
c***********************************************************************
c output parameter:
c
      REAL
     &     velocity(NODE2dVECG(ifirst,ilast,gcw),1:NDIM)
c     
c***********************************************************************
c local parameter
c
      integer i, j
c
c***********************************************************************
c

c      write (6, *) "block_id = ", block_id 
c      do j = ifirst1, ilast1+1
c        do i = ifirst0, ilast0+1
c          write (6, *) "node_coord", i,j, 
c     &                 node_coord(i,j,1), node_coord(i,j,2)
c        enddo
c      enddo

      if (block_id .eq. 0) then
        do j = ifirst1, ilast1+1
          do i = ifirst0, ilast0+1
            velocity(i,j,1) = 0.0d0
            velocity(i,j,2) = 0.0d0
          enddo
        enddo
      else if (block_id .eq. 1) then
        do j = ifirst1, ilast1+1
          do i = ifirst0, ilast0+1
            velocity(i,j,1) = d_velocity_in(1)
            velocity(i,j,2) = d_velocity_in(2)
          enddo
        enddo
      else
        write (6, *) "Wrong block_id = ", block_id
        call flush(6)
      endif

c***********************************************************************
      return
      end   
