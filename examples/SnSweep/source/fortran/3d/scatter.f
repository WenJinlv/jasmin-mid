
c***********************************************************************
      subroutine scatter(
     &  qs,sgm_s,wf,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mg0,mg1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mg0,mg1
      double precision
     &     qs   (0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1),
     &     sgm_s(0:mg1-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1),
     &     wf   (ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,
     &          ifirst2:ilast2+1,0:mg0-1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j,k
      integer mg, ig, lg1
     
      do k = ifirst2, ilast2+1
      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         lg1=0
         do mg=0,mg0-1
            qs(mg,i,j,k)=0d0
            do ig=0,mg
               qs(mg,i,j,k)=qs(mg,i,j,k)+sgm_s(lg1,i,j,k)*wf(i,j,k,ig)
               lg1=lg1+1
            enddo
         enddo
      enddo
      enddo
      enddo

      return
      end   

c***********************************************************************
