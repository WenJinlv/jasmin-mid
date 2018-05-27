define(NDIM,2)dnl
define(REAL,double precision)dnl
include(pdat_m4arrdim2d.i)dnl
c
c     该函数在指定的索引区域,填充网格片影像区.
c
      subroutine physical_Boundary_Conditions(
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  bcs0,bce0,bcs1,bce1,
     &  coord,vel,den,mass,pp,
     &  btype,locindex,bcond)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(const.i)dnl
include(ConstantDefines.h)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer
     &  ics0,ice0,ics1,ice1,icg0,icg1,
     &  bcs0,bce0,bcs1,bce1,
     &  btype,locindex,bcond
c
      REAL  
     &  coord(NODE2dVECG(ics,ice,icg),NDIM),  
     &  vel(NODE2dVECG(ics,ice,icg),NDIM),
     &  den(CELL2dVECG(ics,ice,icg)),
     &  mass(CELL2dVECG(ics,ice,icg)),
     &  pp(CELL2dVECG(ics,ice,icg))
c
      integer 
     &  ist,iend,jst,jend,isym,jsym,isgn,jsgn,sgn,i,j,k,ii,jj,loc
      REAL
     &  p0(NDIM),p1(NDIM),p2(NDIM)
      REAL
     &  a1,r01
c
      if( btype.eq.1.and.
     1    bcond.ne.SYMM_BDY2D_COND_TYPE.and.
     2    bcond.ne.WALL_BDY2D_COND_TYPE.and.
     3    bcond.ne.FREE_BDY2D_COND_TYPE.and.
     4    bcond.ne.DIRI_BDY2D_COND_TYPE.and.
     5    bcond.ne.ZERO_BDY2D_COND_TYPE.and.
     6    bcond.ne.OFLO_BDY2D_COND_TYPE
     7  ) return
c
      if( btype.eq.2.and.
     1    bcond.ne.X_SYMM_CORNER_TYPE.and.
     2    bcond.ne.R_SYMM_CORNER_TYPE.and.
     3    bcond.ne.X_WALL_CORNER_TYPE.and.
     4    bcond.ne.R_WALL_CORNER_TYPE.and.
     5    bcond.ne.X_FREE_CORNER_TYPE.and.
     6    bcond.ne.R_FREE_CORNER_TYPE
     8  ) return
c
      loc=locindex
      if (btype.eq.2) then 
         if( bcond.eq.X_SYMM_CORNER_TYPE.or.
     1       bcond.eq.X_WALL_CORNER_TYPE.or.
     2       bcond.eq.X_FREE_CORNER_TYPE) then
             if(mod(locindex,NDIM).eq.0) then
                loc=0
             else
                loc=1
             endif
         endif
         if( bcond.eq.R_SYMM_CORNER_TYPE.or.
     1       bcond.eq.R_WALL_CORNER_TYPE.or.
     2       bcond.eq.R_FREE_CORNER_TYPE) then
             if(locindex/NDIM.eq.0) then
                loc=2
             else
                loc=3
             endif
         endif
      endif
c
      !填充网格中心量(密度和压力)
      sgn =1
      if(mod(loc,NDIM).eq.1) sgn=-1
      isym=0
      isgn=0
      if(loc.eq.0) isym=bce0
      if(loc.eq.1) isym=bcs0
      if(loc.eq.0.or.loc.eq.1) isgn=1
      jsym=0
      jsgn=0
      if(loc.eq.2) jsym=bce1
      if(loc.eq.3) jsym=bcs1
      if(loc.eq.2.or.loc.eq.3) jsgn=1
      do j=bcs1,bce1
         jj=j
         if(jsgn.eq.1) jj=2*jsym-j+sgn
         do i=bcs0,bce0
            ii=i
            if(isgn.eq.1) ii=2*isym-i+sgn
            den(i,j)=den(ii,jj)
            mass(i,j)=mass(ii,jj)
            pp(i,j)=pp(ii,jj)
            if(bcond.eq.FREE_BDY2D_COND_TYPE.or.
     1         bcond.eq.DIRI_BDY2D_COND_TYPE) then
               den(i,j)=0.0d0
               mass(i,j)=0.0d0
               pp(i,j)=0.0d0
            endif
         enddo
      enddo
c
      !填充网格结点量(坐标和速度)
      isym=0
      jsym=0
      isgn=0
      jsgn=0
      ist =bcs0
      iend=bce0
      jst =bcs1
      jend=bce1
      goto (10,20,30,40) loc+1
 10   isym=bce0+1
      jend=jend+1
      isgn=1
      do k=1,NDIM
         p0(k)=coord(isym,jst ,k)
         p1(k)=coord(isym,jend,k)
      ENDDO
      goto 50
 20   isym=bcs0
      iend=iend+1
      ist =ist+1
      jend=jend+1
      isgn=1
      do k=1,NDIM
         p0(k)=coord(isym,jst ,k)
         p1(k)=coord(isym,jend,k)
      ENDDO
      goto 50
 30   jsym=bce1+1
      iend=iend+1
      jsgn=1
      do k=1,NDIM
         p0(k)=coord(ist ,jsym,k)
         p1(k)=coord(iend,jsym,k)
      ENDDO
      goto 50
 40   jsym=bcs1
      jend=jend+1
      jst=jst+1
      iend=iend+1
      jsgn=1
      do k=1,NDIM
         p0(k)=coord(ist ,jsym,k)
         p1(k)=coord(iend,jsym,k)
      ENDDO
 50   continue
c
      r01=0.0d0
      do k=1,NDIM
         p1(k)=p1(k)-p0(k)
         r01=r01+p1(k)*p1(k)
      enddo
c
      do 100 j=jst,jend
      do 100 i=ist,iend
         ii=i
         if(isgn.eq.1) ii=2*isym-i
         jj=j
         if(jsgn.eq.1) jj=2*jsym-j
         if(bcond.eq.FREE_BDY2D_COND_TYPE.or.
     1      bcond.eq.DIRI_BDY2D_COND_TYPE.or.
     2      bcond.eq.X_FREE_CORNER_TYPE.or.
     3      bcond.eq.R_FREE_CORNER_TYPE) then
c
            if(isgn.eq.1) then
               do k=1,NDIM
                  coord(i,j,k)=2*coord(isym,j,k)-coord(ii,j,k)
c                 vel(i,j,k)  =2*vel(isym,j,k)-vel(ii,j,k)
                  vel(i,j,k)  = 0d0
               enddo
            else if(jsgn.eq.1) then
               do k=1,NDIM
                  coord(i,j,k)=2*coord(i,jsym,k)-coord(i,jj,k)
c                 vel(i,j,k)  =2*vel(i,jsym,k)-vel(i,jj,k)
                  vel(i,j,k)  = 0d0
               enddo
            endif
c
         else
c
            do k=1,NDIM
               p2(k)=coord(ii,jj,k)-p0(k)
            ENDDO
            if(r01<epsilon) then
               a1=0.0d0
            else
               a1=0.0d0
               do k=1,NDIM
                  a1=a1+p1(k)*p2(k)
               enddo
               a1=a1/r01
            endif
            do k=1,NDIM
               coord(i,j,k)=p0(k)+2d0*a1*p1(k)-p2(k)
            enddo
c
            if(r01<epsilon) then
               a1=0.0d0
	    else
               a1=0.0d0
               do k=1,NDIM
                  a1=a1+vel(ii,jj,k)*p1(k)
               enddo
               a1=a1/r01
	    endif
            do k=1,NDIM
               vel(i,j,k)=2d0*a1*p1(k)-vel(ii,jj,k)
            enddo
         endif

 100  continue
c
      return
      end 
