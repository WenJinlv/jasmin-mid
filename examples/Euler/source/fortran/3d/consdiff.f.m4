define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl



c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff(ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  flux0,flux1,flux2,
     &  gamma,density,velocity,pressure)
c***********************************************************************
      implicit none
include(FORTDIR/../probparams.i)dnl
include(FORTDIR/../const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1,ifirst2,ilast2
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     gamma,
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL3d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,ic2,k
      REAL temp,v2norm,mom(NDIM),energy,
     &     gam_min_one
      
c***********************************************************************
c update conserved to full time
c note the permutation of indices in 2nd, 3rd coordinate directions
c***********************************************************************
      gam_min_one = gamma - one

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
           do ic0=ifirst0,ilast0
             mom(1) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,1)
             mom(2) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,2)
             mom(3) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,3)
             v2norm = (velocity(ic0,ic1,ic2,1)**2+
     &                 velocity(ic0,ic1,ic2,2)**2+
     &                 velocity(ic0,ic1,ic2,3)**2)
             energy = pressure(ic0,ic1,ic2)/gam_min_one +
     &                        half*density(ic0,ic1,ic2)*v2norm
 
             density(ic0,ic1,ic2) = density(ic0,ic1,ic2)
     &          -(flux0(ic0+1,ic1,ic2,1)-flux0(ic0,ic1,ic2,1))/dx(0)
     &          -(flux1(ic1+1,ic2,ic0,1)-flux1(ic1,ic2,ic0,1))/dx(1)
     &          -(flux2(ic2+1,ic0,ic1,1)-flux2(ic2,ic0,ic1,1))/dx(2)
             density(ic0,ic1,ic2) = max(smallr,density(ic0,ic1,ic2))

             do k=1,3
               mom(k) = mom(k)
     &          -(flux0(ic0+1,ic1,ic2,k+1)-flux0(ic0,ic1,ic2,k+1))/dx(0)
     &          -(flux1(ic1+1,ic2,ic0,k+1)-flux1(ic1,ic2,ic0,k+1))/dx(1)
     &          -(flux2(ic2+1,ic0,ic1,k+1)-flux2(ic2,ic0,ic1,k+1))/dx(2)
               velocity(ic0,ic1,ic2,k) = mom(k)/density(ic0,ic1,ic2)
             enddo
             energy = energy
     &        -(flux0(ic0+1,ic1,ic2,NEQU)-flux0(ic0,ic1,ic2,NEQU))/dx(0)
     &        -(flux1(ic1+1,ic2,ic0,NEQU)-flux1(ic1,ic2,ic0,NEQU))/dx(1)
     &        -(flux2(ic2+1,ic0,ic1,NEQU)-flux2(ic2,ic0,ic1,NEQU))/dx(2)
c      
             v2norm = (velocity(ic0,ic1,ic2,1)**2+
     &            velocity(ic0,ic1,ic2,2)**2+velocity(ic0,ic1,ic2,3)**2)
             temp = energy - half*density(ic0,ic1,ic2)*v2norm
             pressure(ic0,ic1,ic2) = gam_min_one*temp
             pressure(ic0,ic1,ic2) = max(smallr,pressure(ic0,ic1,ic2))
           enddo
         enddo
      enddo
c
      return
      end
      
  
      subroutine cjconsdiff(model,xlo,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  flux0,flux1,flux2,density,velocity,pressure,df,
     &  gamma,djd,cjd,djv,cje,cjdt,cjnb,lp,rp,time)
c***********************************************************************
      implicit none
include(FORTDIR/../probparams.i)dnl
include(FORTDIR/../const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1, ifirst2, ilast2, model
      REAL xlo(0:NDIM-1),dx(0:NDIM-1),
     &     gamma,djd,cjd,djv,cje,cjdt,cjnb,time,
     &     lp(0:NDIM-1),rp(0:NDIM-1)
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL3d(ifirst,ilast,CELLG)),
     &     df(CELL3d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,ic2,k
      REAL temp,v2norm,mom(NDIM),energy,lr0,rr0,tb,
     &     gam_min_one,f1,f2,dfold,dfnew,dfmax
      REAL xc(0:NDIM-1)
      
c***********************************************************************
c update conserved variables to full time
c***********************************************************************
      gam_min_one = gamma-one 

      do ic2=ifirst2,ilast2
      xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
c
c          璁＄time跺荤у芥板.
           if(model.eq.0) then  ! two-point-detonation
              lr0=sqrt((xc(0)-lp(0))*(xc(0)-lp(0))+
     &                 (xc(1)-lp(1))*(xc(1)-lp(1))+
     &                 (xc(2)-lp(2))*(xc(2)-lp(2)))
              rr0=sqrt((xc(0)-rp(0))*(xc(0)-rp(0))+
     &                 (xc(1)-rp(1))*(xc(1)-rp(1))+
     &                 (xc(2)-rp(2))*(xc(2)-rp(2)))
              tb=lr0/djv
              tb=min(tb,rr0/djv)
           else if(model.eq.1) then  ! plate-detonation
              tb=abs(xc(2)-0.5d0*(lp(2)+rp(2)))/djv
           endif
           if(time.ge.tb+cjdt) then
              dfnew=1.0d0
           else if(time.le.tb) then
              dfnew=0.0d0
           else
              f2=(time-tb)/cjdt
              f1=(1.0d0/djd-1.0d0/density(ic0,ic1,ic2))/
     &           (1.0d0/djd-1.0d0/cjd)
              dfnew=min(1.0d0,(max(f1,f2))**cjnb)
           endif
c
c          璁剧疆璧风′欢, 灏借间负濮瀛.
           dfold=df(ic0,ic1,ic2)
           dfmax=max(dfold,dfnew)
c
           if(dfold.le.1.d-6.and.dfnew.gt.1.d-6) then
              pressure(ic0,ic1,ic2) = gam_min_one*density(ic0,ic1,ic2)
     &                                *cje*dfnew
           else if(dfold.gt.1.d-6) then
              mom(1) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,1)
              mom(2) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,2)
              mom(3) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,3)
              v2norm = ( velocity(ic0,ic1,ic2,1)**2
     &                  +velocity(ic0,ic1,ic2,2)**2 
     &                  +velocity(ic0,ic1,ic2,3)**2)
              energy = pressure(ic0,ic1,ic2)/(gam_min_one*dfold)
     &                       + half*density(ic0,ic1,ic2)*v2norm
 
              density(ic0,ic1,ic2) = density(ic0,ic1,ic2)
     &              -(flux0(ic0+1,ic1,ic2,1)-flux0(ic0,ic1,ic2,1))/dx(0)
     &              -(flux1(ic1+1,ic2,ic0,1)-flux1(ic1,ic2,ic0,1))/dx(1)
     &              -(flux2(ic2+1,ic0,ic1,1)-flux2(ic2,ic0,ic1,1))/dx(2)
              density(ic0,ic1,ic2) = max(smallr,density(ic0,ic1,ic2))
              do k=1,NDIM
                 mom(k) = mom(k)
     &         -(flux0(ic0+1,ic1,ic2,k+1)-flux0(ic0,ic1,ic2,k+1))/dx(0)
     &         -(flux1(ic1+1,ic2,ic0,k+1)-flux1(ic1,ic2,ic0,k+1))/dx(1)
     &         -(flux2(ic2+1,ic0,ic1,k+1)-flux2(ic2,ic0,ic1,k+1))/dx(2)
                 velocity(ic0,ic1,ic2,k) = mom(k)/density(ic0,ic1,ic2)
              enddo
              energy = energy
     &       -(flux0(ic0+1,ic1,ic2,NEQU)-flux0(ic0,ic1,ic2,NEQU))/dx(0)
     &       -(flux1(ic1+1,ic2,ic0,NEQU)-flux1(ic1,ic2,ic0,NEQU))/dx(1)
     &       -(flux2(ic2+1,ic0,ic1,NEQU)-flux2(ic2,ic0,ic1,NEQU))/dx(2)
c      
               v2norm = ( velocity(ic0,ic1,ic2,1)**2
     &                   +velocity(ic0,ic1,ic2,2)**2
     &                   +velocity(ic0,ic1,ic2,3)**2 )
               temp = energy - half*density(ic0,ic1,ic2)*v2norm
               pressure(ic0,ic1,ic2) = gam_min_one*temp*dfmax
               pressure(ic0,ic1,ic2) = max(smallr,pressure(ic0,ic1,ic2))
           endif
c
           df(ic0,ic1,ic2) = dfmax
c
c      if(xc(2).gt.dx(2)) then
c      if(pressure(ic0,ic1,ic2).ge.37.2.and.
c    &    density(ic0,ic1,ic2).ge.2.48) then
c             print *,"Touch the maximum at time=",time,
c    &    " pres,density,eng,vel,sound,df=",pressure(ic0,ic1,ic2),
c    &    density(ic0,ic1,ic2),temp/density(ic0,ic1,ic2),
c    &    velocity(ic0,ic1,ic2,2),
c    &    sqrt(gamma*pressure(ic0,ic1,ic2)/density(ic0,ic1,ic2)),
c    &               df(ic0,ic1,ic2)
c             stop
c          endif
c       endif
c
        enddo
      enddo
      enddo
c
      return
      end      
c***********************************************************************

