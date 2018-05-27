define(NDIM,2)dnl
define(NEQU,1)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine fluxcorrec(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dx,
     &  uval,
     &  flux0,flux1,
     &  trlft0,trlft1,
     &  trrgt0,trrgt1)
c***********************************************************************
      implicit none
include(FORTDIR/../const.i)dnl
include(FORTDIR/../probparams.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      REAL dt 
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     uval(CELL2d(ifirst,ilast,CELLG)),
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)), 
     &     trlft0(FACE2d0(ifirst,ilast,FACEG)),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG)),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG)),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL trnsvers
c    
      do ic1=ifirst1-1,ilast1+1
        do ic0=ifirst0-1,ilast0+1
          trnsvers= (flux0(ic0+1,ic1)-flux0(ic0,ic1))*0.5/dx(0)

          trrgt1(ic1  ,ic0)= trrgt1(ic1  ,ic0) - trnsvers
          trlft1(ic1+1,ic0)= trlft1(ic1+1,ic0) - trnsvers
        enddo
      enddo
c
      do ic0=ifirst0-1,ilast0+1
        do ic1=ifirst1-1,ilast1+1
          trnsvers= (flux1(ic1+1,ic0)-flux1(ic1,ic0))*0.5/dx(1)
          trrgt0(ic0  ,ic1)= trrgt0(ic0  ,ic1) - trnsvers
          trlft0(ic0+1,ic1)= trlft0(ic0+1,ic1) - trnsvers
        enddo
      enddo
c
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcalculation(dt,extra_cell,visco,dx,xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  uval,
     &  flux0,flux1,
     &  trlft0,trlft1,trrgt0,trrgt1)
c***********************************************************************
      implicit none
include(FORTDIR/../const.i)dnl
include(FORTDIR/../probparams.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,extra_cell,visco
      REAL dt 
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     uval(CELL2d(ifirst,ilast,CELLG)),
c variables in 2d side indexed         
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)), 
     &     trlft0(FACE2d0(ifirst,ilast,FACEG)),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG)),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG)),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1,ie0,ie1
      REAL   riemst,vel
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c
      do ic1=ifirst1,ilast1
      vel=-(xlo(1)+dx(1)*(ic1-ifirst1+half))
      do ie0=ifirst0,ilast0+1
         if (vel.ge.zero) then
            riemst= trlft0(ie0,ic1)
         else
            riemst= trrgt0(ie0,ic1)
         endif
         flux0(ie0,ic1)= dt*riemst*vel
      enddo
      enddo      

      do ic0=ifirst0,ilast0
      vel=xlo(0)+dx(0)*(ic0-ifirst0+half)
      do ie1=ifirst1,ilast1+1
         if (vel.ge.zero) then
               riemst= trlft1(ie1,ic0)
         else
               riemst= trrgt1(ie1,ic0)
         endif
         flux1(ie1,ic0)= dt*riemst*vel
      enddo
      enddo
c
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff(ifirst0,ilast0,ifirst1,ilast1,
     &  dx,
     &  flux0,flux1,
     &  uval)
c***********************************************************************
      implicit none
include(FORTDIR/../const.i)dnl
include(FORTDIR/../probparams.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)),
     &     uval(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1
c     
c***********************************************************************
      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0
          uval(ic0,ic1) = uval(ic0,ic1)
     &          -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &          -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
        enddo
      enddo
c
      return
      end
