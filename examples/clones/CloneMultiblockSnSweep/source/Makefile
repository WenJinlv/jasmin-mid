##
## 文件名:	Makefile.in
## 软件包:	JASMIN applications
## 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
## 版本号:	$Revision: 1.6 $
## 修改  :	$Date: 2007/05/21 01:22:16 $
## 描述  :	makefile for rigid rotation sample application
##

JASMIN	      = /home/wenjinlv/jasmin-mid

FORTRAN       = fortran
M4DIRS        = -I $(FORTRAN)/$(NDIM)d -I $(JASMIN)/include

default:      main-2d

include $(JASMIN)/lib/jasmin/config/Makefile.config

CPPFLAGS_EXTRA= -DNDIM=$(NDIM) 
#CPPFLAGS_EXTRA= -DNDIM=$(NDIM) -DTESTING

CXX_OBJS      = SnSweep.o  SnSweepLevelIntegrator.o  Main.o 
FOR_OBJS      = init.o     scatter.o  solve.o   sourcew.o  vflux.o   computeEnergy.o

main-2d:
		if test -f stamp-3d; then $(MAKE) clean ; fi
		touch stamp-2d
		$(MAKE) NDIM=2 main2d

main2d:		$(CXX_OBJS) $(FOR_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) \
		$(CXX_OBJS) $(FOR_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main2d

main-3d:
		if test -f stamp-2d; then $(MAKE) clean ; fi
		touch stamp-3d
		$(MAKE) NDIM=3 main3d

main3d:		$(CXX_OBJS) $(FOR_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) \
		$(CXX_OBJS) $(FOR_OBJS)	\
		$(LIBJASMIN3D) $(LIBJASMIN) $(CXXLDLIBS) -o main3d
clean: 
		$(JASCLEAN)
		$(RM) *.f main main2d main3d


