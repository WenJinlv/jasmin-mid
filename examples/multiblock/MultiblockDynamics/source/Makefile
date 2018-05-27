##
## 文件名:	Makefile.in
## 软件包:	JASMIN applications
## 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
## 版本号:	$Revision: 170 $
## 修改  :	$Date: 2007-06-27 08:44:22  $
## 描述  :	makefile for Euler gas dynamics sample application
##

JASMIN	      = /home/wenjinlv/jasmin-mid

FORTRAN       = fortran
M4DIRS        = -I fortran/$(NDIM)d -I $(JASMIN)/include

          

default:        main-2d

include $(JASMIN)/lib/jasmin/config/Makefile.config

CPPFLAGS_EXTRA= -DNDIM=$(NDIM) -DERROR_DIFFUSE=$(ERROR_DIFFUSE) -DERROR_ROUNDOFF=$(ERROR_ROUNDOFF)

CXX_OBJS      = Dynamics.o  DynamicsLevelIntegrator.o  Main.o
F2D_OBJS      = alternalNodeCoupler.o       computeSpeedup.o      \
                physicalBoundaryConditions.o  stableDt.o          \
                computeGridDensityEnergy.o  computeVelocity.o     postprocessVelocity.o \
                computeSidePressure.o       initializeVelocity.o  preprocessDynamicsState.o \
                pside.o negvol.o posvol.o 


main-2d: main-error-roundoff
	mv main_error_roundoff main2d


main-error-diffuse:
		$(MAKE) clean
		$(MAKE) NDIM=2 ERROR_DIFFUSE=1 ERROR_ROUNDOFF=0 main_error_diffuse

main_error_diffuse:		$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) $(CXX_OBJS) $(F2D_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main_error_diffuse

main-error-roundoff:
		$(MAKE) clean
		$(MAKE) NDIM=2 ERROR_DIFFUSE=1 ERROR_ROUNDOFF=1 main_error_roundoff

main_error_roundoff:		$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) $(CXX_OBJS) $(F2D_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main_error_roundoff

main-no-error-diffuse:
		$(MAKE) clean
		$(MAKE) NDIM=2 ERROR_DIFFUSE=0 ERROR_ROUNDOFF=0 main_no_error_diffuse

main_no_error_diffuse:		$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) $(CXX_OBJS) $(F2D_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main_no_error_diffuse

clean:
		$(RM) *.F *.o stamp-* 

redo:
		$(RM) core *.ii *.int.c 
		$(RM) main2d main_no_error_diffuse main_error_roundoff main_error_diffuse
		$(RM) -r ti_files ii_files 


