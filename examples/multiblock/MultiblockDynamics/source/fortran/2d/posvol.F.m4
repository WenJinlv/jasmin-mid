
define(REAL,double precision)dnl

      function posvol(x1,r1,x2,r2,x3,r3,x4,r4)
      implicit none

      REAL 
     1 x1,r1,x2,r2,x3,r3,x4,r4,posvol

include(const.i)dnl
include(ConstantDefines.h)dnl
c                                                                      
      posvol= THIRD*(x1-x2)*(r1*r1+r1*r2+r2*r2)
     1    + THIRD*(x2-x3)*(r2*r2+r2*r3+r3*r3)
     2    + THIRD*(x3-x4)*(r3*r3+r3*r4+r4*r4)
     3    + THIRD*(x4-x1)*(r4*r4+r4*r1+r1*r1)

      if(posvol.lt.ZERO.and.r1.le.0.0d0.and.r2.le.0.0d0.and.
     1   r3.le.0.0d0.and.r4.le.0.0d0) posvol=-posvol
	                                                      
      return
      end
