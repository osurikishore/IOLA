include ../../configure.hwrf

EXEC = ../../exec/hwrf_da_update_analysis.exe

MODS = kinds_interface.o \
       constants_interface.o \
       namelist_interface.o \
       diagnostics_interface.o \
       variable_interface.o \
       grid_methods_interface.o \
       math_methods_interface.o \
       time_methods_interface.o \
       netcdf_interface.o \
       json_interface.o \
       fileio_interface.o \
       interpolation_interface.o \
       bndrycnd_interface.o \
       nmm_interpolation_interface.o \
       analysis_interface.o \
       hwrf_analysis_interface.o \
       da_update_analysis_interface.o

SUBS =

OBJS = $(MODS) $(SUBS) main.o

INC = -I../../libs/mods/fson -I../../libs/mods/slint -I../../libs/mods/kdtree $(NETCDF_INCLUDE) -I$(MKLROOT)/include/fftw
LIBS = -L../../libs -lspline -lslatec -lfson -lslint -lkdtree $(NETCDF_LDFLAGS)
FFLAGS += $(INC) -mkl=sequential -DLINUX

$(EXEC): main.F90 $(MODS) $(SUBS)
	$(SFC) $(FFLAGS) main.F90 $(MODS) $(SUBS) -o $(@) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(MODS) *.mod *.f90

.F90.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.F90
