












c***********************************************************************
c
c    Initialization routine 
c
c***********************************************************************
      subroutine init(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mg0, mg1, ms0,
     &  sgm_t, sgm_s, 
     &  wf, fi, 
     &  dx, xleng,
     &  sgm_t0, sgm_s0,
     &  vn,
     &  angles, sw,
     &  am, bm)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mg0, mg1, ms0
      double precision
     &     sgm_t(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1),
     &     sgm_s(0:mg1-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1),
     &     wf   (ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:mg0-1),
     &     fi(0:mg0-1,ifirst0:ilast0+1,
     &          ifirst1:ilast1+1,0:ms0-1)
      double precision 
     &     dx(0:2-1),xleng(0:2-1),
     &     sgm_t0(0:mg0-1), sgm_s0(0:mg1-1),
     &     vn(0:mg0-1),
     &     angles(0:2-1,0:ms0-1),sw(0:ms0-1)
      double precision 
     &     am(ifirst0:ilast0+1,0:1),
     &     bm(ifirst1:ilast1+1,0:1)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *     
      integer i,j, ms,mg,ig,lg1, l0,l1

c     *** initialize sgm_s and sgm_t ***

      do j = ifirst1, ilast1+1
      do i = ifirst0, ilast0+1
         lg1 = 0
         do mg= 0, mg0-1
            sgm_t(mg,i,j)=sgm_t0(mg)
            do ig=0,mg
               sgm_s(lg1,i,j)=sgm_s0(lg1) 
               lg1=lg1+1
            enddo
         enddo
      enddo
      enddo

c     *** initialize fi ***

      
      do i = ifirst0, ilast0+1
         am(i,0) = i*dx(0)/xleng(0)
         am(i,1) = 1d0 - am(i,0)
      enddo
      do j = ifirst1, ilast1+1
         bm(j,0) = j*dx(1)/xleng(1)
         bm(j,1) = 1d0 - bm(j,0)
      enddo


      do ms = 0, ms0-1
         l0 = 0
         l1 = 0
         if ( angles(0,ms) .lt. 0.d0 ) l0 = 1
         if ( angles(1,ms) .lt. 0.d0 ) l1 = 1
         do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0+1
         do mg= 0, mg0-1
            fi(mg,i,j,ms) = vn(mg) * am(i,l0)*bm(j,l1)
         enddo
         enddo
         enddo
      enddo

c     *** initialize wf ***
      call vflux(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     mg0, ms0,
     &     wf, fi, 
     &     sw)

      return
      end   

c***********************************************************************

      subroutine snquad(angles,sw,u0,ms0,n0)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

      implicit none
      integer ms0, n0
      double precision
     &     angles(0:2-1,0:ms0-1),sw(0:ms0-1),
     &     u0(0:n0/2-1)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

      integer n2, ms4, ms, i, j

      n2=n0/2        
      ms4=ms0/4

      ms=0
      do i=0,n2-1
         do j=0,n2-i-1
            angles(0,ms)=u0(i)
            angles(1,ms)=u0(j)
            sw(ms)=sw(ms)
            angles(0,ms+ms4)=-u0(i)
            angles(1,ms+ms4)=u0(j)
            sw(ms+ms4)=sw(ms)
            angles(0,ms+2*ms4)=-u0(i)
            angles(1,ms+2*ms4)=-u0(j)
            sw(ms+2*ms4)=sw(ms)
            angles(0,ms+3*ms4)=u0(i)
            angles(1,ms+3*ms4)=-u0(j)
            sw(ms+3*ms4)=sw(ms)
            ms=ms+1
         enddo
      enddo

      return
      end

c***********************************************************************
