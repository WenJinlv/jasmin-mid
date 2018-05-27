
include(JASMIN_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine tagpatchforgeom(xlo_global, xup_global, dx,xlo,
     &                    ifirst0,ilast0,
     &                    ifirst1,ilast1,
     &                    ifirst2,ilast2,
     &                    gwd0,gwd1,gwd2,
     &                    tag)
      implicit none
      real*8 dx(0:2), xlo(0:2), xlo_global(0:2), xup_global(0:2)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        ifirst2,ilast2,
     &        gwd0,gwd1, gwd2

      integer tag(CELL3d(ifirst,ilast,0))

      integer i,j,k
      real*8  loc(0:2),xwid(0:2)

      xwid(0)=xup_global(0)-xlo_global(0)
      xwid(1)=xup_global(1)-xlo_global(1)
      xwid(2)=xup_global(2)-xlo_global(2)
      do k=ifirst2,ilast2
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         loc(0)=xlo(0)+(i-ifirst0+0.5d0)*dx(0)-xlo_global(0)
         loc(1)=xlo(1)+(j-ifirst1+0.5d0)*dx(1)-xlo_global(1)
         loc(2)=xlo(2)+(k-ifirst1+0.5d0)*dx(2)-xlo_global(2)
         tag(i,j,k)=0
         if(loc(0).lt.xwid(0)*0.75d0) tag(i,j,k)=1
      enddo
      enddo
      enddo

      return
      end
