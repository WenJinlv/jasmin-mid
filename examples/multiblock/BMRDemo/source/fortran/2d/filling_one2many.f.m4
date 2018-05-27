define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    填充一对多影像区: 边型.
c
c***********************************************************************
      subroutine filling_one2many(
     &  ifirst_g0,ilast_g0,ifirst_g1,ilast_g1,
     &  gcwg0, gcwg1,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0, gcw1,
     &  ifirst_f0,ilast_f0,ifirst_f1,ilast_f1,
     &  ratio_in0, ratio_in1,
     &  center_var_ghost, center_var)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst_f0,ilast_f0,ifirst_f1,ilast_f1
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ifirst_g0,ilast_g0,ifirst_g1,ilast_g1
      integer gcw0, gcw1, gcwg0, gcwg1
      integer ratio_in0, ratio_in1
c
c variables in 2d indexed
c
      REAL
     &     center_var(CELL2dVECG(ifirst,ilast,gcw)),
     &     center_var_ghost(CELL2dVECG(ifirst_g,ilast_g,gcwg))
c
      REAL center_var_local
      integer ratio0, ratio1
      integer i,j, i_g,j_g, i_g_beg, j_g_beg, i_d,j_d, i_r,j_r
      logical is_print
c
c***********************************************************************     
c

      
      if (ratio_in0 >1) then

      ratio0 = ratio_in0
      ratio1 = ratio_in1
      j_d = ifirst_g1*ratio0-ifirst_f1
      i_d = ifirst_g0*ratio0-ifirst_f0
      do j = ifirst_f1, ilast_f1
         j_g = (j+j_d)/ratio1
         do i = ifirst_f0, ilast_f0 
            i_g = (i+i_d)/ratio0
            center_var(i,j)=center_var_ghost(i_g,j_g)
         enddo
      enddo
 
      else if (ratio_in0 <1) then

      ratio0 = -ratio_in0
      ratio1 = -ratio_in1
c      print*,"ratio_in0, ratio_in1= ", ratio_in0, ratio_in1

      j_d = ifirst_g1-ifirst_f1*ratio1
      i_d = ifirst_g0-ifirst_f0*ratio0

      do j = ifirst_f1, ilast_f1
      j_g_beg = j*ratio1
      do i = ifirst_f0, ilast_f0
         i_g_beg = i*ratio0
         center_var_local = 0.0d0

c         print*,"i,j= ", i,j
c         print*,"ratio0, ratio1= ", ratio0,ratio1

         do j_r = 0, ratio1-1
         j_g = j_g_beg+j_r+j_d
         do i_r = 0, ratio0-1
            i_g = i_g_beg+i_r+i_d

c            print*, "i_g, j_g= ", i_g, j_g
c            print*, "center_var_ghost= ", center_var_ghost(i_g,j_g)

            center_var_local=center_var_local+center_var_ghost(i_g,j_g)
         enddo
         enddo

c         print*,"center_var_local= ", center_var_local

         center_var_local = center_var_local/(ratio0*ratio1)
         center_var(i,j) = center_var_local

c         print*,"center_var(i,j)= ", center_var(i,j)
      enddo
      enddo
 

      else 

      j_d = ifirst_g1-ifirst_f1
      i_d = ifirst_g0-ifirst_f0
      do j = ifirst_f1, ilast_f1
         j_g = j+j_d
         do i = ifirst_f0, ilast_f0 
            i_g = i+i_d
            center_var(i,j)=center_var_ghost(i_g,j_g)
         enddo
      enddo
 
      endif

c
      return
      end

c***********************************************************************
c
c    填充一对多影像区: 结点型.
c
c***********************************************************************
      subroutine filling_one2many_node(
     &  ifirst_g0,ilast_g0,ifirst_g1,ilast_g1,
     &  gcwg0, gcwg1,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0, gcw1,
     &  ifirst_f0,ilast_f0,ifirst_f1,ilast_f1,
     &  ratio_in0, ratio_in1,
     &  center_var_ghost, center_var)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst_f0,ilast_f0,ifirst_f1,ilast_f1
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ifirst_g0,ilast_g0,ifirst_g1,ilast_g1
      integer gcw0, gcw1, gcwg0, gcwg1
      integer ratio_in0, ratio_in1
c
c variables in 2d indexed
c
      REAL
     &     center_var(CELL2dVECG(ifirst,ilast,gcw)),
     &     center_var_ghost(CELL2dVECG(ifirst_g,ilast_g,gcwg))
c
      REAL center_var_local
      integer ratio0, ratio1
      integer i,j, i_g,j_g, i_g_beg, j_g_beg, i_d,j_d, i_r,j_r
      logical is_print
c
c***********************************************************************     
c

      
      if (ratio_in0 >1) then

      ratio0 = ratio_in0
      ratio1 = ratio_in1
      j_d = ifirst_g1*ratio0-ifirst_f1
      i_d = ifirst_g0*ratio0-ifirst_f0

      ilast_f0 = ifirst_f0+gcw0-1
      ilast_f1 = ifirst_f1+gcw1-1

      do j = ifirst_f1, ilast_f1
         j_g = (j+j_d)/ratio1
         do i = ifirst_f0, ilast_f0 
            i_g = (i+i_d)/ratio0
            center_var(i,j)=center_var(i,j)+center_var_ghost(i_g,j_g)
         enddo
      enddo
 
      else if (ratio_in0 <1) then

      ratio0 = -ratio_in0
      ratio1 = -ratio_in1
c      print*,"ratio_in0, ratio_in1= ", ratio_in0, ratio_in1

      j_d = ifirst_g1-ifirst_f1*ratio1
      i_d = ifirst_g0-ifirst_f0*ratio0

      ilast_f0 = ifirst_f0+gcw0-1
      ilast_f1 = ifirst_f1+gcw1-1

      do j = ifirst_f1, ilast_f1
      j_g_beg = j*ratio1
      do i = ifirst_f0, ilast_f0
         i_g_beg = i*ratio0
         center_var_local = 0.0d0

c         print*,"i,j= ", i,j
c         print*,"ratio0, ratio1= ", ratio0,ratio1

         do j_r = 0, ratio1-1
         j_g = j_g_beg+j_r+j_d
         do i_r = 0, ratio0-1
            i_g = i_g_beg+i_r+i_d

c            print*, "i_g, j_g= ", i_g, j_g
c            print*, "center_var_ghost= ", center_var_ghost(i_g,j_g)

            center_var_local=center_var_local+center_var_ghost(i_g,j_g)
         enddo
         enddo

c         print*,"center_var_local= ", center_var_local

         center_var_local = center_var_local/(ratio0*ratio1)
         center_var(i,j) = center_var(i,j)+center_var_local

c         print*,"center_var(i,j)= ", center_var(i,j)
      enddo
      enddo
 

      else 

      j_d = ifirst_g1-ifirst_f1
      i_d = ifirst_g0-ifirst_f0
      do j = ifirst_f1, ilast_f1
         j_g = j+j_d
         do i = ifirst_f0, ilast_f0 
            i_g = i+i_d
            center_var(i,j)=center_var(i,j)+center_var_ghost(i_g,j_g)
         enddo
      enddo
 
      endif

c
      return
      end




