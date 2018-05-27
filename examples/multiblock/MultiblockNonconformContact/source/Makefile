##
## 文件名:	Makefile.in
## 软件包:	JASMIN applications
## 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
## 版本号:	$Revision: 1.5 $
## 修改  :	$Date: 2007/08/09 02:01:06 $
## 描述  :	makefile for Dynamic application
##

JASMIN          = /home/wenjinlv/jasmin-mid

FORTRAN       = fortran
M4DIRS        = -DFORTDIR=$(FORTRAN) -DJASMIN_FORTDIR=$(JASMIN)/include

default:        main-2d

include $(JASMIN)/lib/jasmin/config/Makefile.config

CPPFLAGS_EXTRA = -DNDIM=$(NDIM) 
FFLAGS_EXTRA   =
#FFLAGS_EXTRA   = -fno-automatic
#LDFLAGS_EXTRA= --link_command_prefix "/usr/local/bin/purify4.5 -best-effort"

CXX_OBJS = Dynamic.o DynamicLevelIntegrator.o main_Dynamic.o 

F2D_OBJS      = initialize_field_2block_horizontal.o \
                initialize_field_2block_vertical.o \
                initialize_field_3block_horizontal.o \
                advance_coordinate.o 

main-2d:
		if test -f stamp-3d; then $(MAKE) clean; fi
		touch stamp-2d
		$(MAKE) NDIM=2 main2d

main2d:	$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) $(F2D_OBJS) \
		$(LIBJASMIN2D) $(LIBJASMIN) $(LDLIBS) -o main2d 

clean:
		$(RM) *.mod *.f *.f90 *.o main2d main2ds

cleanall: 
	$(RM) -rf viz* restart* log.* stamp-*d
	$(RM) -rf *.mod *.f *.f90 *.o main2d main3d *.orig core*


