include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c     
c  更新结点位置坐标.
c     
c***********************************************************************
      subroutine advancecoordinate(
     &     ifirst0, ilast0,
     &     ifirst1, ilast1,
     &     gcw0, gcw1,
     &     dt,
     &     velocity,
     &     node_coord)
c***********************************************************************
      implicit none
c***********************************************************************
c input parameter:
c
      integer ifirst0, ilast0, ifirst1, ilast1,
     &        gcw0, gcw1

      REAL dt,
     &     velocity(NODE2dVECG(ifirst,ilast,gcw),1:NDIM)
c
c***********************************************************************
c output parameter:
c
      REAL 
     &     node_coord(NODE2dVECG(ifirst,ilast,gcw),1:NDIM)
c     
c***********************************************************************
c local parameter
c
      integer i, j, k
c
c***********************************************************************
c
      do j = ifirst1, ilast1 + 1
        do i = ifirst0, ilast0 + 1
c          write (6,*) i,j,velocity(i,j,1),velocity(i,j,2)
c          write (6,*) i,j,node_coord(i,j,1),node_coord(i,j,2)
c          call flush(6)
          do k = 1, NDIM
            ! 更新结点坐标.
            node_coord(i,j,k) = node_coord(i,j,k) + velocity(i,j,k) * dt
          enddo
c          write (6,*) i,j,node_coord(i,j,1),node_coord(i,j,2)
c          call flush(6)
        enddo
      enddo

c***********************************************************************
      return
      end   
