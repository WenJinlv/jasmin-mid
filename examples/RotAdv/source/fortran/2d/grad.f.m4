define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl
 
      subroutine detectgrad(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  dx,
     &  gradtol,
     &  dotag,donttag,
     &  var,
     &  tags)
c***********************************************************************
      implicit none
include(FORTDIR/../probparams.i)dnl
include(FORTDIR/../const.i)dnl
c***********************************************************************
c input arrays:
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  dotag,donttag,
     &  vghost0,vghost1,
     &  tagghost0,tagghost1,
     &  ttagghost0,ttagghost1
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
c variables indexed as 2dimensional
      REAL
     &  var(CELL2dVECG(ifirst,ilast,vghost))
      integer
     &  tags(CELL2dVECG(ifirst,ilast,tagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1,diag01
      logical tagcell
      integer ic0,ic1
c
c***********************************************************************
c
      tol = gradtol
      diag01 = sqrt(dx(0)**2+dx(1)**2)

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0

          if (tags(ic0,ic1) .ne. 0) then
            loctol = 0.125*tol
          else 
            loctol = tol
          endif
   
          tagcell = .false.

          presm1 = var(ic0-1,ic1)
          presp1 = var(ic0+1,ic1)
          facejump = abs(var(ic0,ic1)-presm1)
          facejump = max(facejump,abs(var(ic0,ic1)-presp1))
          tagcell = ((facejump).gt.(loctol*dx(0)))
          if (.not.tagcell) then
            presm1 = var(ic0,ic1-1)
            presp1 = var(ic0,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*dx(1)))
          endif

          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1-1)
            presp1 = var(ic0+1,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif
          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1+1)
            presp1 = var(ic0+1,ic1-1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif

          if ( tagcell ) then
            tags(ic0,ic1) = dotag
          endif
        enddo
      enddo
c
      return
      end


