c***********************************************************************
c
c   根据物理边界条件，填充影像区数据.
c
c***********************************************************************
      subroutine fillphysicbdry(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  igfirst0,iglast0,igfirst1,iglast1,
     &  bloc, 
     &  gcw0,gcw1,
     &  clone_number,
     &  uval)

      implicit none

c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1     ! 网格片索引范围
      integer igfirst0,iglast0,igfirst1,iglast1 ! 物理边界影像区索引范围
      integer bloc                              ! 物理边界影像区位置编号
      integer gcw0,gcw1,clone_number            ! uval的数据片影像区宽度
      double precision
     &     uval(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1)

c  JASMIN框架中与边界位置相关的常量定义
      include "BoundaryDefines.i" 

c  
      integer i,j
      real*8  ratio

      ratio=clone_number*1.0d0;
      if(clone_number.le.0) ratio=0.0d0 
      if(bloc.eq.XLO) then
          do j=igfirst1,iglast1
          do i=igfirst0,iglast0
               uval(i,j)=2.0d0*(ratio+1)
          enddo
          enddo
      else if(bloc.eq.XHI) then
          do j=igfirst1,iglast1
          do i=igfirst0,iglast0
               uval(i,j)=0.0d0
          enddo
          enddo
      endif
 
      return
      end
