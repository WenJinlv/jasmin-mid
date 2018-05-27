define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    设置边界条件: 网格结点量.
c
c***********************************************************************
      subroutine setphysbdryfornodes(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  xlo,xhi,ylo,yhi,xdx,ydy,
     &  igfirst0,iglast0,igfirst1,iglast1,
     &  type,loc,xo,yo)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,
     &        igfirst0,iglast0,igfirst1,iglast1
      integer gcw0,gcw1,type,loc
      REAL
     &        xlo,xhi,ylo,yhi,xdx,ydy
c
c variables in 2d indexed
c
      REAL
     &     xo(NODE2dVECG(ifirst,ilast,gcw)),
     &     yo(NODE2dVECG(ifirst,ilast,gcw))
c
      integer i,j,l
c
c***********************************************************************     
c
      if(type.eq.1) then
         if(loc.eq.0) then
            do j=igfirst1,iglast1+1
            do i=igfirst0,iglast0
               xo(i,j)=xlo-(i-igfirst0+1)*xdx
               yo(i,j)=ylo+j*ydy
            enddo
            enddo
         else if(loc.eq.1) then
            do j=igfirst1  ,iglast1+1
            do i=igfirst0+1,iglast0+1
               xo(i,j)=xhi+(i-igfirst0)*xdx
               yo(i,j)=ylo+j*ydy
            enddo
            enddo
         else if(loc.eq.2) then
            do j=igfirst1,iglast1
            do i=igfirst0,iglast0+1
               xo(i,j)=xlo+i*xdx
               yo(i,j)=ylo-(j-igfirst1+1)*ydy
            enddo
            enddo
         else 
            do j=igfirst1+1,iglast1+1
            do i=igfirst0  ,iglast0+1
               xo(i,j)=xlo+i*xdx
               yo(i,j)=yhi+(j-igfirst1)*ydy
            enddo
            enddo
         endif
      endif
c
      if(type.eq.2) then
         if(loc.eq.0) then
            xo(igfirst0,igfirst1)=xo(igfirst0+1,igfirst1)-xdx
            yo(igfirst0,igfirst1)=yo(igfirst0+1,igfirst1)
         else if(loc.eq.1) then
            xo(iglast0+1,igfirst1)=xo(iglast0,igfirst1)+xdx
            yo(iglast0+1,igfirst1)=yo(iglast0,igfirst1)
         else if(loc.eq.2) then
            xo(igfirst0,iglast1+1)=xo(igfirst0+1,iglast1+1)-xdx
            yo(igfirst0,iglast1+1)=yo(igfirst0+1,iglast1+1)
         else
            xo(iglast0+1,iglast1+1)=xo(iglast0,iglast1+1)+xdx
            yo(iglast0+1,iglast1+1)=yo(iglast0,iglast1+1)
         endif
      endif
c
      return
      end
