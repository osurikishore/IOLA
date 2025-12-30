include ../../configure.hwrf

EXEC = ../../exec/hwrf_tempdrop.exe

MODS = kinds_interface.o \
       namelist_interface.o \
       constants_interface.o \
       diagnostics_interface.o \
       variable_interface.o \
       netcdf_interface.o \
       fileio_interface.o \
       time_methods_interface.o \
       grid_methods_interface.o \
       math_methods_interface.o \
       meteo_methods_interface.o \
       tempdrop_sonde_interface.o 

SUBS =

OBJS = $(MODS) $(SUBS) main.o

INC = $(NETCDF_INCLUDE)
LIBS = -L../../libs -lslatec -lnr -lsondelib $(NETCDF_LDFLAGS)
FFLAGS += $(INC)

$(EXEC): main.F90 $(MODS) $(SUBS)
	$(SFC) $(FFLAGS) main.F90 $(MODS) $(SUBS) -o $(@) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(MODS) *.mod *.f90

.F90.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.F90