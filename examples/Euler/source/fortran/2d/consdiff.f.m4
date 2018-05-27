define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl


c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff(ifirst0,ilast0,ifirst1,ilast1,dx,
     &  flux0,flux1,
     &  gamma,density,velocity,pressure)
c***********************************************************************
      implicit none
include(FORTDIR/../probparams.i)dnl
include(FORTDIR/../const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG),NEQU),
     &     gamma,
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,k
      REAL temp,v2norm,mom(NDIM),energy,
     &     gam_min_one
      
c***********************************************************************
c update conserved variables to full time
c note the permutation of indices in 2nd coordinate direction
c***********************************************************************
      gam_min_one = gamma-one 

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0
          mom(1) = density(ic0,ic1)*velocity(ic0,ic1,1)
          mom(2) = density(ic0,ic1)*velocity(ic0,ic1,2)
          v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
          energy = pressure(ic0,ic1)/gam_min_one +
     &                        half*density(ic0,ic1)*v2norm
 
          density(ic0,ic1) = density(ic0,ic1)
     &          -(flux0(ic0+1,ic1,1)-flux0(ic0,ic1,1))/dx(0)
     &          -(flux1(ic1+1,ic0,1)-flux1(ic1,ic0,1))/dx(1)
          density(ic0,ic1) = max(smallr,density(ic0,ic1))
          do k=1,NDIM
            mom(k) = mom(k)
     &        -(flux0(ic0+1,ic1,k+1)-flux0(ic0,ic1,k+1))/dx(0)
     &        -(flux1(ic1+1,ic0,k+1)-flux1(ic1,ic0,k+1))/dx(1)
            velocity(ic0,ic1,k) = mom(k)/density(ic0,ic1)
          enddo
          energy = energy
     &       -(flux0(ic0+1,ic1,NEQU)-flux0(ic0,ic1,NEQU))/dx(0)
     &       -(flux1(ic1+1,ic0,NEQU)-flux1(ic1,ic0,NEQU))/dx(1)
c      
          v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
          temp = energy - half*density(ic0,ic1)*v2norm
          pressure(ic0,ic1) = gam_min_one*temp
          pressure(ic0,ic1) = max(smallr,pressure(ic0,ic1))
        enddo
      enddo
c
      return
      end
      
       subroutine cjconsdiff(model,xlo,dx,ifirst0,ilast0,ifirst1,ilast1,
     &  flux0,flux1,density,velocity,pressure,df,
     &  gamma,djd,cjd,djv,cje,cjdt,cjnb,lp,rp,time)
c***********************************************************************
      implicit none
include(FORTDIR/../probparams.i)dnl
include(FORTDIR/../const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1, model
      REAL xlo(0:NDIM-1),dx(0:NDIM-1),
     &     gamma,djd,cjd,djv,cje,cjdt,cjnb,time,
     &     lp(0:NDIM-1),rp(0:NDIM-1)
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG),NEQU),
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL2d(ifirst,ilast,CELLG)),
     &     df(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,k
      REAL temp,v2norm,mom(NDIM),energy,lr0,rr0,tb,
     &     gam_min_one,f1,f2,dfold,dfnew,dfmax
      REAL xc(0:NDIM-1)
      
c***********************************************************************
c update conserved variables to full time
c note the permutation of indices in 2nd coordinate direction
c***********************************************************************
      gam_min_one = gamma-one 

      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
c
c          璁＄time跺荤у芥板.
           if(model.eq.0) then  ! two-point-detonation
              lr0=sqrt((xc(0)-lp(0))*(xc(0)-lp(0))+
     &                 (xc(1)-lp(1))*(xc(1)-lp(1)))
              rr0=sqrt((xc(0)-rp(0))*(xc(0)-rp(0))+
     &                 (xc(1)-rp(1))*(xc(1)-rp(1)))
              tb=lr0/djv
              tb=min(tb,rr0/djv)
           else if(model.eq.1) then  ! plate-detonation
              tb=abs(xc(1)-0.5d0*(lp(1)+rp(1)))/djv
           endif
           if(time.ge.tb+cjdt) then
              dfnew=1.0d0
           else if(time.le.tb) then
              dfnew=0.0d0
           else
              f2=(time-tb)/cjdt
              f1=(1.0d0/djd-1.0d0/density(ic0,ic1))/
     &           (1.0d0/djd-1.0d0/cjd)
              dfnew=min(1.0d0,(max(f1,f2))**cjnb)
           endif
c
c          璁剧疆璧风′欢, 灏借间负濮瀛.
           dfold=df(ic0,ic1)
           dfmax=max(dfold,dfnew)
c
           if(dfold.le.1.d-6.and.dfnew.gt.1.d-6) then
              pressure(ic0,ic1) = gam_min_one*density(ic0,ic1)
     &                              *cje*dfnew
           else if(dfold.gt.1.d-6) then
              mom(1) = density(ic0,ic1)*velocity(ic0,ic1,1)
              mom(2) = density(ic0,ic1)*velocity(ic0,ic1,2)
              v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
              energy = pressure(ic0,ic1)/(gam_min_one*dfold)
     &                       + half*density(ic0,ic1)*v2norm
 
              density(ic0,ic1) = density(ic0,ic1)
     &              -(flux0(ic0+1,ic1,1)-flux0(ic0,ic1,1))/dx(0)
     &              -(flux1(ic1+1,ic0,1)-flux1(ic1,ic0,1))/dx(1)
              density(ic0,ic1) = max(smallr,density(ic0,ic1))
              do k=1,NDIM
                 mom(k) = mom(k)
     &              -(flux0(ic0+1,ic1,k+1)-flux0(ic0,ic1,k+1))/dx(0)
     &              -(flux1(ic1+1,ic0,k+1)-flux1(ic1,ic0,k+1))/dx(1)
                 velocity(ic0,ic1,k) = mom(k)/density(ic0,ic1)
              enddo
              energy = energy
     &             -(flux0(ic0+1,ic1,NEQU)-flux0(ic0,ic1,NEQU))/dx(0)
     &             -(flux1(ic1+1,ic0,NEQU)-flux1(ic1,ic0,NEQU))/dx(1)
c      
               v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
               temp = energy - half*density(ic0,ic1)*v2norm
               pressure(ic0,ic1) = gam_min_one*temp*dfmax
               pressure(ic0,ic1) = max(smallr,pressure(ic0,ic1))
           endif
c
           df(ic0,ic1) = dfmax
c
c      if(xc(1).gt.dx(1)) then
c      if(pressure(ic0,ic1).ge.37.2.and.density(ic0,ic1).ge.2.48) then
c             print *,"Touch the maximum at time=",time,
c    &    " pres,density,eng,vel,sound,df=",pressure(ic0,ic1),
c    &    density(ic0,ic1),temp/density(ic0,ic1),
c    &    velocity(ic0,ic1,2),
c    &    sqrt(gamma*pressure(ic0,ic1)/density(ic0,ic1)),df(ic0,ic1)
c             stop
c          endif
c       endif
c
        enddo
      enddo
c
      return
      end
c***********************************************************************

