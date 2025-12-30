include ../../configure.hwrf

EXEC = ../../exec/hwrf_obstobufr.exe

MODS = kinds_interface.o \
       namelist_interface.o \
       constants_interface.o \
       diagnostics_interface.o \
       variable_interface.o \
       bufrio_interface.o \
       fileio_interface.o \
       time_methods_interface.o \
       math_methods_interface.o \
       meteo_methods_interface.o \
       obs_to_bufr_interface.o 

SUBS =

OBJS = $(MODS) $(SUBS) main.o

INC = $(NETCDF_INCLUDE)
LIBS = -L../../libs -lslatec -lnr -lbufr_i4r8
FFLAGS += $(INC)

$(EXEC): main.F90 $(MODS) $(SUBS)
	$(SFC) $(FFLAGS) main.F90 $(MODS) $(SUBS) -o $(@) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(MODS) *.mod *.f90

.F90.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.F90