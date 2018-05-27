define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    设置物理边界条件: 网格单元中心量.
c
c***********************************************************************
      subroutine setphysbdryforcells(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  igfirst0,iglast0,igfirst1,iglast1,
     &  type,loc,lmax,uo)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,
     &        igfirst0,iglast0,igfirst1,iglast1
      integer gcw0,gcw1,type,loc,lmax
c
c variables in 2d indexed
c
      REAL
     &     uo(CELL2dVECG(ifirst,ilast,gcw),lmax)
c
      integer i,j,l
c
c***********************************************************************     
C
      if(type.eq.1) then
         if(loc.eq.0) then
            do j=igfirst1,iglast1
            do i=igfirst0,iglast0
               do l=1,lmax
                  uo(i,j,l)=uo(iglast0+1+(iglast0-i),j,l)
               enddo
            enddo
            enddo
         else if(loc.eq.1) then
            do j=igfirst1,iglast1
            do i=igfirst0,iglast0
               do l=1,lmax
                  uo(i,j,l)=uo(igfirst0-1-(i-igfirst0),j,l)
               enddo
            enddo
            enddo
         else if(loc.eq.2) then
            do j=igfirst1,iglast1
            do i=igfirst0,iglast0
               do l=1,lmax
                  uo(i,j,l)=uo(i,iglast1+1+(iglast1-j),l)
               enddo
            enddo
            enddo
         else 
            do j=igfirst1,iglast1
            do i=igfirst0,iglast0
               do l=1,lmax
                  uo(i,j,l)=uo(i,igfirst1-1-(j-igfirst1),l)
               enddo
            enddo
            enddo
         endif
      endif
c
      if(type.eq.2) then
         if(loc.eq.0) then
            do l=1,lmax
               uo(igfirst0,igfirst1,l)=uo(igfirst0+1,igfirst1,l)
            enddo
         else if(loc.eq.1) then
            do l=1,lmax
               uo(iglast0,igfirst1,l)=uo(iglast0-1,igfirst1,l)
            enddo
         else if(loc.eq.2) then
            do l=1,lmax
               uo(igfirst0,iglast1,l)=uo(igfirst0+1,iglast1,l)
            enddo
         else
            do l=1,lmax
               uo(iglast0,iglast1,l)=uo(iglast0-1,iglast1,l)
            enddo
         endif
      endif
c
      return
      end
