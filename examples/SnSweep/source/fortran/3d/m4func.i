define(CAL_AM_BM_CM,`dnl

      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo
      do k = ifirst2, ilast2+1
         cm(k,0) = k*dx(2)/xleng(2)
         cm(k,1) = 1d0 - cm(k,0)
      enddo
')dnl

