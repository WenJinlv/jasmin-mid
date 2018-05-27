define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    BMR插值: 不同类型的变量采用不同的子程序.
c
c***********************************************************************

c=======================================================================
c    BMR插值: 坐标.
c=======================================================================
      subroutine interpolating_coord( ifirst_src0, ilast_src0,
     &                          ifirst_src1, ilast_src1,
     &                          ghost_bmr0, ghost_bmr1,
     &                          ifirst_dst0, ilast_dst0,
     &                          ifirst_dst1, ilast_dst1,
     &                          ratio_dst_to_src0, ratio_dst_to_src1,
     &                          coords_bmr,
     &                          coords_current)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst_src0, ilast_src0, ifirst_src1, ilast_src1
      integer ifirst_dst0, ilast_dst0, ifirst_dst1, ilast_dst1
      integer ghost_bmr0, ghost_bmr1
      integer ratio_dst_to_src0, ratio_dst_to_src1
c
c variables in 2d indexed
c
      REAL
     & coords_current(NODE2d(ifirst_dst,ilast_dst,0),0:NDIM-1),
     & coords_bmr(NODE2dVECG(ifirst_src,ilast_src,ghost_bmr),0:NDIM-1)
c
      integer i,j,k, ii,jj,iii, k0,k1
      REAL ratio, dist0, dist1, ratio_dtos0, ratio_dtos1
      integer ki, kj
      REAL ratio0, ratio1
c
c***********************************************************************     
c
      REAL, allocatable :: coords_temp(:,:,:)
      integer ibeg0, iend0, ibeg1, iend1, error

      ratio_dtos0 = ratio_dst_to_src0
      ratio_dtos1 = ratio_dst_to_src1

c========================================
c     细化插值.
c========================================
      if (ratio_dst_to_src0 .gt. 1) then

c        create temporal variable
c================================
         ibeg0 = ifirst_src0*ratio_dst_to_src0
         iend0 = (ilast_src0+1)*ratio_dst_to_src0-1
         ibeg1 = ifirst_src1*ratio_dst_to_src1
         iend1 = (ilast_src1+1)*ratio_dst_to_src1-1
         allocate( coords_temp(ibeg0:iend0+1,ibeg1:iend1+1,0:NDIM-1),
     &             stat = error)
         call chkallocate(error)

c        overlap points
c======================
         DO j=ifirst_src1,ilast_src1+1
         DO i=ifirst_src0,ilast_src0+1
            jj = j*ratio_dst_to_src1
            ii = i*ratio_dst_to_src0
            coords_temp(ii,jj,0) = coords_bmr(i,j,0)
            coords_temp(ii,jj,1) = coords_bmr(i,j,1)
         ENDDO
         ENDDO
c        new points on old x-line
c=================================
         DO j=ifirst_src1,ilast_src1+1
         DO i=ifirst_src0,ilast_src0
            dist0 = coords_bmr(i+1,j,0)-coords_bmr(i,j,0)
            dist1 = coords_bmr(i+1,j,1)-coords_bmr(i,j,1)
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k0 = 1,ratio_dst_to_src0-1
               ratio = k0/ratio_dtos0
               coords_temp(ii+k0,jj,0)=coords_temp(ii,jj,0)+
     &                                    dist0*ratio
               coords_temp(ii+k0,jj,1)=coords_temp(ii,jj,1)+
     &                                    dist1*ratio
            enddo
         ENDDO
         ENDDO
c        new points on old y-line
c=================================
         DO i=ifirst_src0,ilast_src0+1
         DO j=ifirst_src1,ilast_src1
            dist0 = coords_bmr(i,j+1,0)-coords_bmr(i,j,0)
            dist1 = coords_bmr(i,j+1,1)-coords_bmr(i,j,1)
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k1 = 1,ratio_dst_to_src1-1
               ratio = k1/ratio_dtos1
               coords_temp(ii,jj+k1,0)=coords_temp(ii,jj,0)+
     &                                    dist0*ratio
               coords_temp(ii,jj+k1,1)=coords_temp(ii,jj,1)+
     &                                    dist1*ratio
            enddo
         ENDDO
         ENDDO
c        new points on new lines
c=================================
         DO j=ifirst_src1,ilast_src1
         DO i=ifirst_src0,ilast_src0
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k1 = 1,ratio_dst_to_src1-1
               iii = ii+ratio_dst_to_src0
               dist0 = coords_temp(iii,jj+k1,0)-
     &                 coords_temp(ii,jj+k1,0)
               dist1 = coords_temp(iii,jj+k1,1)-
     &                 coords_temp(ii,jj+k1,1)
               do k0 = 1,ratio_dst_to_src0-1
                  ratio = k0/ratio_dtos0
                  coords_temp(ii+k0,jj+k1,0)=
     &               coords_temp(ii,jj+k1,0)+dist0*ratio
                  coords_temp(ii+k0,jj+k1,1)=
     &               coords_temp(ii,jj+k1,1)+dist1*ratio
               enddo
            enddo
         ENDDO
         ENDDO

         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            coords_current(i,j,0) = coords_temp(i,j,0)
            coords_current(i,j,1) = coords_temp(i,j,1)
         ENDDO
         ENDDO

         deallocate(coords_temp)
c========================================
c     粗化插值.
c========================================
      else if (ratio_dst_to_src0 .lt. 0) then
         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            jj = -j*ratio_dst_to_src1
            ii = -i*ratio_dst_to_src0
            coords_current(i,j,0) = coords_bmr(ii,jj,0)
            coords_current(i,j,1) = coords_bmr(ii,jj,1)
         ENDDO
         ENDDO

c========================================
c     直接复制.
c========================================
      else 
         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            coords_current(i,j,0) = coords_bmr(i,j,0)
            coords_current(i,j,1) = coords_bmr(i,j,1)
         ENDDO
         ENDDO
      endif

c     check
      DO j=ifirst_dst1,ilast_dst1+1
      DO i=ifirst_dst0,ilast_dst0+1
         if (coords_current(i,j,0).gt.100000.0) then
            print*, "i,j = ", i,j
            print*, "coords_current(i,j,0)=", coords_current(i,j,0)
         endif
         if (coords_current(i,j,1).gt.100000.0) then
            print*, "i,j = ", i,j
            print*, "coords_current(i,j,1)=", coords_current(i,j,1)
         endif
      ENDDO
      ENDDO
c

      return
      end


c=======================================================================
c    BMR插值: 中心型.
c=======================================================================
      subroutine interpolating_center( ifirst_src0, ilast_src0,
     &                          ifirst_src1, ilast_src1,
     &                          ghost_bmr0, ghost_bmr1,
     &                          ifirst_dst0, ilast_dst0,
     &                          ifirst_dst1, ilast_dst1,
     &                          ratio_dst_to_src0, ratio_dst_to_src1,
     &                          source_bmr,
     &                          source,
     &                          center_var_bmr,
     &                          center_var_current)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst_src0, ilast_src0, ifirst_src1, ilast_src1
      integer ifirst_dst0, ilast_dst0, ifirst_dst1, ilast_dst1
      integer ghost_bmr0, ghost_bmr1
      integer ratio_dst_to_src0, ratio_dst_to_src1
c
c variables in 2d indexed
c
      REAL
     &     source(CELL2d(ifirst_dst,ilast_dst,0)),
     &     source_bmr(CELL2d(ifirst_src,ilast_src,0)),
     &     center_var_current(CELL2d(ifirst_dst,ilast_dst,0)),
     &     center_var_bmr(CELL2dVECG(ifirst_src,ilast_src,ghost_bmr))
c
      integer i,j,ii,jj,k0,k1
      logical is_print
c
c***********************************************************************     
c
      integer ibeg0, iend0, ibeg1, iend1, error
      REAL, allocatable :: center_var_temp(:,:)

      is_print = .false.

c========================================
c     细化插值.
c========================================
      if (ratio_dst_to_src0 .gt. 1) then

c        create temporal variable
c================================
         ibeg0 = ifirst_src0*ratio_dst_to_src0
         iend0 = (ilast_src0+1)*ratio_dst_to_src0-1
         ibeg1 = ifirst_src1*ratio_dst_to_src1
         iend1 = (ilast_src1+1)*ratio_dst_to_src1-1
         allocate( center_var_temp(ibeg0:iend0,ibeg1:iend1),
     &             stat = error)
         call chkallocate(error)

c        intepolation
c================================
         DO j=ifirst_src1,ilast_src1
         DO i=ifirst_src0,ilast_src0

         if (is_print) then
               if (center_var_bmr(i-1,j).gt.1.0e5) then
                  print*, "i-1,j= ", i-1, j, ", center_var_bmr=", 
     &                    center_var_bmr(i-1,j)
               endif
               if (center_var_bmr(i+1,j).gt.1.0e5) then
                  print*, "i+1,j= ", i+1, j, ", center_var_bmr=", 
     &                    center_var_bmr(i+1,j)
               endif
               if (center_var_bmr(i,j-1).gt.1.0e5) then
                  print*, "i,j-1= ", i, j-1, ", center_var_bmr=", 
     &                    center_var_bmr(i,j-1)
               endif
               if (center_var_bmr(i,j+1).gt.1.0e5) then
                  print*, "i,j+1= ", i, j+1, ", center_var_bmr=", 
     &                    center_var_bmr(i,j+1)
               endif
         endif

            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k1 = 0,ratio_dst_to_src1-1
            do k0 = 0,ratio_dst_to_src0-1

               center_var_temp(ii+k0,jj+k1)=1.0*center_var_bmr(i,j)
c     &                                       +0.1*center_var_bmr(i-1,j)
c     &                                       +0.1*center_var_bmr(i+1,j)
c     &                                       +0.1*center_var_bmr(i,j-1)
c     &                                       +0.1*center_var_bmr(i,j+1)
            enddo
            enddo
         ENDDO
         ENDDO

         DO j=ifirst_dst1,ilast_dst1
         DO i=ifirst_dst0,ilast_dst0
            center_var_current(i,j) = center_var_temp(i,j)
         ENDDO
         ENDDO

         deallocate(center_var_temp)

c========================================
c     粗化插值.
c========================================
      else if (ratio_dst_to_src0 .lt. 0) then
         DO j=ifirst_dst1,ilast_dst1
         DO i=ifirst_dst0,ilast_dst0
            jj = -j*ratio_dst_to_src1
            ii = -i*ratio_dst_to_src0
            center_var_current(i,j) = center_var_bmr(ii,jj)
         ENDDO
         ENDDO

c========================================
c     直接复制.
c========================================
      else 
         DO j=ifirst_dst1,ilast_dst1
         DO i=ifirst_dst0,ilast_dst0
            center_var_current(i,j) = center_var_bmr(i,j)
         ENDDO
         ENDDO
      endif

c
      return
      end


c=======================================================================
c    BMR插值: 结点型.
c=======================================================================
      subroutine interpolating_node( ifirst_src0, ilast_src0,
     &                          ifirst_src1, ilast_src1,
     &                          ghost_bmr0, ghost_bmr1,
     &                          ifirst_dst0, ilast_dst0,
     &                          ifirst_dst1, ilast_dst1,
     &                          ratio_dst_to_src0, ratio_dst_to_src1,
     &                          node_var_bmr,
     &                          node_var_current)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst_src0, ilast_src0, ifirst_src1, ilast_src1
      integer ifirst_dst0, ilast_dst0, ifirst_dst1, ilast_dst1
      integer ghost_bmr0, ghost_bmr1
      integer ratio_dst_to_src0, ratio_dst_to_src1
c
c variables in 2d indexed
c
      REAL
     &     node_var_current(NODE2d(ifirst_dst,ilast_dst,0)),
     &     node_var_bmr(NODE2dVECG(ifirst_src,ilast_src,ghost_bmr))
c
      integer i,j,ii,jj,iii,jjj,k0,k1
      REAL ratiox_l, ratiox_r, ratioy_l, ratioy_r
      REAL var_l, var_r, var_ll, var_lr, var_rl, var_rr
      REAL vol_ll, vol_lr, vol_rl, vol_rr
      REAL ratio_dtos0, ratio_dtos1
c
c***********************************************************************     
c
      integer ibeg0, iend0, ibeg1, iend1, error
      REAL, allocatable :: node_var_temp(:,:)

      ratio_dtos0 = ratio_dst_to_src0
      ratio_dtos1 = ratio_dst_to_src1

c========================================
c     细化插值.
c========================================
      if (ratio_dst_to_src0 .gt. 1) then

c        create temporal variable
c================================
         ibeg0 = ifirst_src0*ratio_dst_to_src0
         iend0 = (ilast_src0+1)*ratio_dst_to_src0-1
         ibeg1 = ifirst_src1*ratio_dst_to_src1
         iend1 = (ilast_src1+1)*ratio_dst_to_src1-1
         allocate( node_var_temp(ibeg0:iend0+1,ibeg1:iend1+1),
     &             stat = error)
         call chkallocate(error)

c        overlap points
c======================
         DO j=ifirst_src1,ilast_src1+1
         DO i=ifirst_src0,ilast_src0+1
            jj = j*ratio_dst_to_src1
            ii = i*ratio_dst_to_src0
            node_var_temp(ii,jj) = node_var_bmr(i,j)
         ENDDO
         ENDDO
c        new points on old x-line
c=================================
         DO j=ifirst_src1,ilast_src1+1
         DO i=ifirst_src0,ilast_src0
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k0 = 1,ratio_dst_to_src0-1
               ratiox_l = k0/ratio_dtos0
               ratiox_r = 1.0-ratiox_l
               var_l = node_var_temp(ii,jj)
               var_r = node_var_temp(ii+ratio_dst_to_src0,jj)
               node_var_temp(ii+k0,jj)=ratiox_l*var_l+
     &                                    ratiox_r*var_r
            enddo
         ENDDO
         ENDDO
c        new points on old y-line
c=================================
         DO i=ifirst_src0,ilast_src0+1
         DO j=ifirst_src1,ilast_src1
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            do k1 = 1,ratio_dst_to_src1-1
               ratioy_l = k1/ratio_dtos1
               ratioy_r = 1.0-ratioy_l
               var_l = node_var_temp(ii,jj)
               var_r = node_var_temp(ii,jj+ratio_dst_to_src1)
               node_var_temp(ii,jj+k1)=ratioy_l*var_l+
     &                                    ratioy_r*var_r
            enddo
         ENDDO
         ENDDO
c        new points on new lines
c=================================
         DO j=ifirst_src1,ilast_src1
         DO i=ifirst_src0,ilast_src0
            ii = i*ratio_dst_to_src0
            jj = j*ratio_dst_to_src1
            iii = (i+1)*ratio_dst_to_src0
            jjj = (j+1)*ratio_dst_to_src1
            var_ll = node_var_temp(ii,jj)
            var_lr = node_var_temp(ii,jjj)
            var_rl = node_var_temp(iii,jj)
            var_rr = node_var_temp(iii,jjj)
            do k1 = 1,ratio_dst_to_src1-1
               ratioy_l = k1/ratio_dtos1
               ratioy_r = 1.0-ratioy_l
               do k0 = 1,ratio_dst_to_src0-1
                  ratiox_l = k0/ratio_dtos0
                  ratiox_r = 1.0-ratiox_l
                  vol_ll = ratiox_l*ratioy_l
                  vol_lr = ratiox_l*ratioy_r           
                  vol_rl = ratiox_r*ratioy_l
                  vol_rr = ratiox_r*ratioy_r
                  node_var_temp(ii+k0,jj+k1)=
     &                            vol_ll*var_ll+vol_lr*var_lr+
     &                            vol_rl*var_rl+vol_rr*var_rr
               enddo
            enddo
         ENDDO
         ENDDO

         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            node_var_current(i,j) = node_var_temp(i,j)
         ENDDO
         ENDDO

         deallocate( node_var_temp )

c========================================
c     粗化插值.
c========================================
      else if (ratio_dst_to_src0 .lt. 0) then
         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            jj = -j*ratio_dst_to_src1
            ii = -i*ratio_dst_to_src0
            node_var_current(i,j) = node_var_bmr(ii,jj)
         ENDDO
         ENDDO

c========================================
c     直接复制.
c========================================
      else 
         DO j=ifirst_dst1,ilast_dst1+1
         DO i=ifirst_dst0,ilast_dst0+1
            node_var_current(i,j) = node_var_bmr(i,j)
         ENDDO
         ENDDO
      endif
c

      return
      end


c***********************************************************************
c     
c  检查内存开辟函数的返回码是否正确.
c     
c***********************************************************************
      subroutine chkallocate(error)
c***********************************************************************
      implicit none
c***********************************************************************
c input parameter:
c
      integer error
c***********************************************************************
      if (error .ne. 0) then
        write (6, *) "allocate error"
        stop
      endif
c***********************************************************************
      return
      end   

