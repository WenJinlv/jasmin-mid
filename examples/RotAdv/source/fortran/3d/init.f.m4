define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

c***********************************************************************
c
c    Initialization routine
c
c***********************************************************************
      subroutine initrotation(model,dx,xlo,center,radius,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  uval)
c***********************************************************************
      implicit none
include(FORTDIR/../const.i)dnl
include(FORTDIR/../probparams.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer model
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),
     &     center(0:NDIM-1),radius
c variables in 2d cell indexed         
      REAL
     &     uval(CELL3dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer d,ic0,ic1,ic2
      REAL 
     &xc(0:NDIM-1),scenter(0:NDIM-1),r0,r1,r2,sr0,sr1,sr2,
     &angle0,angle1,angle2

c     write(6,*) "in initrotation"
c     write(6,*) "model_problem= ",model
c     write(6,*) "dx = ",(dx(i),i=0,NDIM-1)
c     write(6,*) "xlo = ",(xlo(i),i=0,NDIM-1)
c     write(6,*) "ce = ",(center(i),i=0,NDIM-1)
c     write(6,*) "radius = ",radius
c     write(6,*) "ifirst0,ilast0 = ",ifirst0,ilast0
c     write(6,*) "ifirst1,ilast1 = ",ifirst1,ilast1
c     write(6,*) "ifirst2,ilast2 = ",ifirst2,ilast2
c
      do d=0,NDIM-1
         scenter(d)=-center(d)
      enddo
      do ic2=ifirst2,ilast2
         xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
         r2 = xc(2)-center(2)
         sr2= xc(2)-scenter(2)
         do ic1=ifirst1,ilast1
            xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
            r1 = xc(1)-center(1)
            sr1= xc(1)-scenter(1)
            do ic0=ifirst0,ilast0
               xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
               r0 = xc(0)-center(0)
               sr0= xc(0)-scenter(0)
               if(model.eq.CIRCLE_ROTATION_DELTA) then
                  if(sqrt(r0*r0+r1*r1+r2*r2).le.radius.or.
     &               sqrt(sr0*sr0+sr1*sr1+sr2*sr2).le.radius) then
                     uval(ic0,ic1,ic2)=1.0d0
                  else
                     uval(ic0,ic1,ic2)=0.0d0
                  endif
               else if(model.eq.CIRCLE_ROTATION_PARAB) then
                  angle0 = r0*r0+1.5*r1*r1+1.5*r2*r2
                  angle1 = sr0*sr0+1.5*sr1*sr1+1.5*sr2*sr2
                  if(angle0.le.1.0d0/16) then
                     uval(ic0,ic1,ic2)=1.0d0-16*angle0
                  else if(angle1.le.1.0d0/16) then
                     uval(ic0,ic1,ic2)=1.0d0-16*angle1
                  else
                     uval(ic0,ic1,ic2)=0.0d0
                  endif
               else if(model.eq.SQUARE_ROTATION_DELTA) then
                  if((abs(r0).le.radius.and.abs(r1).le.radius)
     &               .and.abs(r2).le.radius.or.
     &               (abs(sr0).le.radius.and.abs(sr1).le.radius)
     &               .and.abs(sr2).le.radius) then
                     uval(ic0,ic1,ic2)=1.0d0
                  else
                     uval(ic0,ic1,ic2)=0.0d0
                  endif
               endif
            enddo
         enddo
      enddo
c
      return
      end
