define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    Initialization routine where we use a spherical profile 
c
c***********************************************************************
      subroutine initrotation(model,dx,xlo,center,radius,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  uval)
c***********************************************************************
      implicit none
include(FORTDIR/../const.i)dnl
include(FORTDIR/../probparams.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer model
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),
     &     center(0:NDIM-1),radius
c variables in 2d cell indexed         
      REAL
     &     uval(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer d,ic0,ic1
      REAL 
     &xc(0:NDIM-1),scenter(0:NDIM-1),r0,r1,sr0,sr1,angle0,angle1

c     write(6,*) "in initrotation"
c     write(6,*) "model_problem= ",model
c     write(6,*) "dx = ",(dx(i),i=0,NDIM-1)
c     write(6,*) "xlo = ",(xlo(i),i=0,NDIM-1)
c     write(6,*) "ce = ",(center(i),i=0,NDIM-1)
c     write(6,*) "radius = ",radius
c     write(6,*) "ifirst0,ilast0 = ",ifirst0,ilast0
c     write(6,*) "ifirst1,ilast1 = ",ifirst1,ilast1
c
      do d=0,NDIM-1
         scenter(d)=-center(d)
      enddo
      do ic1=ifirst1,ilast1
         xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
         r1 = xc(1)-center(1)
         sr1= xc(1)-scenter(1)
         do ic0=ifirst0,ilast0
            xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
            r0 = xc(0)-center(0)
            sr0= xc(0)-scenter(0)
            if(model.eq.CIRCLE_ROTATION_DELTA) then
               if(sqrt(r0*r0+r1*r1).le.radius.or.
     &            sqrt(sr0*sr0+sr1*sr1).le.radius) then
                  uval(ic0,ic1)=1.0d0
               else
                  uval(ic0,ic1)=0.0d0
               endif
            else if(model.eq.CIRCLE_ROTATION_PARAB) then
               angle0 = r0*r0+1.5*r1*r1
               angle1 = sr0*sr0+1.5*sr1*sr1
               if(angle0.le.1.0d0/16) then
                  uval(ic0,ic1)=1.0d0-16*angle0
               else if(angle1.le.1.0d0/16) then
                  uval(ic0,ic1)=1.0d0-16*angle1
               else
                  uval(ic0,ic1)=0.0d0
               endif
            else if(model.eq.SQUARE_ROTATION_DELTA) then
               if((abs(r0).le.radius.and.abs(r1).le.radius).or.
     &            (abs(sr0).le.radius.and.abs(sr1).le.radius)) then
                  uval(ic0,ic1)=1.0d0
               else
                  uval(ic0,ic1)=0.0d0
               endif
            endif
         enddo
      enddo
c
      return
      end
