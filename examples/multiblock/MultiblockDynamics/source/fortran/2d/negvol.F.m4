
define(REAL,double precision)dnl

      function negvol(x1,r1,x2,r2,x3,r3,x4,r4)
      implicit none

      REAL 
     1 x1,r1,x2,r2,x3,r3,x4,r4,negvol

include(const.i)dnl
include(ConstantDefines.h)dnl
c                                                                      
      negvol= 0.5d0*(1-LCY)*(x1-x2)*(r1+r2) + 
     1           THIRD*LCY*(x1-x2)*(r1*r1+r1*r2+r2*r2)
     2    +0.5d0*(1-LCY)*(x2-x3)*(r2+r3) +
     1           THIRD*LCY*(x2-x3)*(r2*r2+r2*r3+r3*r3)
     2    +0.5d0*(1-LCY)*(x3-x4)*(r3+r4) +
     1           THIRD*LCY*(x3-x4)*(r3*r3+r3*r4+r4*r4)
     2    +0.5d0*(1-LCY)*(x4-x1)*(r4+r1) +
     1           THIRD*LCY*(x4-x1)*(r4*r4+r4*r1+r1*r1)

      return
      end
