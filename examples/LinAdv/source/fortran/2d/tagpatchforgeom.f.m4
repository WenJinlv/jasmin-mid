
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine tagpatchforgeom(xlo_global, xup_global, dx,xlo,
     &                    ifirst0,ilast0,
     &                    ifirst1,ilast1,
     &                    gwd0,gwd1,
     &                    tag)
      implicit none
      real*8 dx(0:1), xlo(0:1), xlo_global(0:1), xup_global(0:1)
      integer ifirst0,ilast0,
     &        ifirst1,ilast1,
     &        gwd0,gwd1

      integer tag(CELL2d(ifirst,ilast,0))

      integer i,j
      real*8  loc(0:1),xwid(0:1)

      xwid(0)=xup_global(0)-xlo_global(0)
      xwid(1)=xup_global(1)-xlo_global(1)
      do j=ifirst1,ilast1
      do i=ifirst0,ilast0
         loc(0)=xlo(0)+(i-ifirst0+0.5d0)*dx(0)-xlo_global(0)
         loc(1)=xlo(1)+(j-ifirst1+0.5d0)*dx(1)-xlo_global(1)
         tag(i,j)=0
         if(loc(0).lt.xwid(0)*0.75d0) tag(i,j)=1
      enddo
      enddo

      return
      end
