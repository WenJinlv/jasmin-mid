c***********************************************************************
c
c   根据物理边界条件，填充影像区数据.
c
c***********************************************************************
      subroutine fillphysicbdry(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  igfirst0,iglast0,igfirst1,iglast1,igfirst2,iglast2,
     &  bloc, 
     &  gcw0,gcw1,gcw2,
     &  clone_number,
     &  uval)

      implicit none

c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2  ! 网格片索引范围
      integer igfirst0,iglast0,igfirst1,iglast1,igfirst2,iglast2 ! 边界影像区
      integer bloc                              ! 物理边界影像区位置编号
      integer gcw0,gcw1,gcw2,clone_number       ! uval的数据片影像区宽度
      double precision
     &     uval(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,
     &          ifirst2-gcw2:ilast2+gcw2)

c  JASMIN框架中与边界位置相关的常量定义
      include "BoundaryDefines.i" 

c  
      integer i,j,k
      real*8  ratio
 
      ratio=clone_number*1.0d0;
      if(clone_number.le.0) ratio=0.0d0 
      if(bloc.eq.XLO) then
          do k=igfirst2,iglast2
          do j=igfirst1,iglast1
          do i=igfirst0,iglast0
               uval(i,j,k)=2.0d0*(ratio+1)
          enddo
          enddo
          enddo
      else if(bloc.eq.XHI) then
          do k=igfirst2,iglast2
          do j=igfirst1,iglast1
          do i=igfirst0,iglast0
               uval(i,j,k)=0.0d0
          enddo
          enddo
          enddo
      endif
 
      return
      end
