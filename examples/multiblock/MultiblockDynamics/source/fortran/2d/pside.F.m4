
define(NDIM,2)dnl
define(REAL,double precision)dnl

      subroutine pside(den,pres,x,r,p0)
      implicit none
c
      REAL 
     1 den,pres,x(2*NDIM),r(2*NDIM),p0(2*NDIM),
     2 ddd,area,f1,f2,f3,f4,a4,a40,bb0,
     3 p1,p2,p3,p4,x1,x2,x3,x4,r1,r2,r3,r4,
     4 a1,b1,a2,b2,a3,b3
c
      ddd(a1,b1,a2,b2,a3,b3)=
     1           0.5d0*(a1*b2+a2*b3+a3*b1-a1*b3-a2*b1-a3*b2)
c                                                                      
      x1=x(1)
      x2=x(2)
      x3=x(3)
      x4=x(4)
      r1=r(1)
      r2=r(2)
      r3=r(3)
      r4=r(4)
c
      area=((x3-x1)*(r4-r2)+(x2-x4)*(r3-r1))*0.5d0
      f1=ddd(x4,r4,x1,r1,x2,r2)
      f2=ddd(x1,r1,x2,r2,x3,r3)
      f3=ddd(x2,r2,x3,r3,x4,r4)
      f4=ddd(x3,r3,x4,r4,x1,r1)
      a4=area/4d0
      a40=area/40d0
      bb0=a4/a40
      if(f1.gt.a4) p1=pres
      if(f1.lt.a4.and.f1.ge.a40) p1=pres*a4/f1
      if(f1.lt.a40) p1=pres*bb0

      if(f2.gt.a4) p2=pres
      if(f2.lt.a4.and.f2.ge.a40) p2=pres*a4/f2
      if(f2.lt.a40) p2=pres*bb0

      if(f3.gt.a4) p3=pres
      if(f3.lt.a4.and.f3.ge.a40) p3=pres*a4/f3
      if(f3.lt.a40) p3=pres*bb0

      if(f4.gt.a4) p4=pres
      if(f4.lt.a4.and.f4.ge.a40) p4=pres*a4/f4
      if(f4.lt.a40) p4=pres*bb0

      p0(1)=0.5d0*(p1+p2)
      p0(2)=0.5d0*(p2+p3)
      p0(3)=0.5d0*(p3+p4)
      p0(4)=0.5d0*(p4+p1)
c
      return
      end
