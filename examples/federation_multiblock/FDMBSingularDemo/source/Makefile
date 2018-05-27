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
M4DIRS        = -DFORTDIR=$(FORTRAN)/$(NDIM)d -DJASMIN_FORTDIR=$(JASMIN)/include

default:      main-2d

include $(JASMIN)/lib/jasmin/config/Makefile.config

CPPFLAGS_EXTRA= -DNDIM=$(NDIM)

CXX_OBJS      = main_FDMBSingDemo.o  FDMBSingDemo.o \
                FDMBSingDemoLevelIntegrator.o \
                FDMBSingDemoFederationLevelIntegrator.o \
                FDMBSingDemoFederation.o
F2D_OBJS      = summing.o

main-2d:
		if test -f stamp-3d; then $(MAKE) clean ; fi
		touch stamp-2d
		$(MAKE) NDIM=2 main2d

main2d:		$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) $(CXX_OBJS) $(F2D_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main2d

clean: 
		$(JASCLEAN)
		$(RM) *.f main main2d main3d


