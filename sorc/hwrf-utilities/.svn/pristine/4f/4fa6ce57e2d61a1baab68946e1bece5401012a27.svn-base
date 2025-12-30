include ../../configure.hwrf

EXEC = ../../exec/hwrf_ensemble.exe

# Define modules and associated libraries

MODS = kinds.o                                       \
       constants.o                                   \
       mpi-interface.o                               \
       namelist.o                                    \
       common-interface.o                            \
       kdtree2.o                                     \
       interpolation-interface.o                     \
       ensemble-mean-interface.o                     \
       ensemble-variance-interface.o                 \
       ensemble-recenter-interface.o                 \
       ensemble-inflate-interface.o                  \
       ensemble-update-interface.o                   \
       wrf-ensemble-interface.o
SUBS =
OBJS=$(MODS) $(SUBS) Driver.o

INC = $(NETCDF_INCLUDE)
LIBS = $(NETCDF_LDFLAGS)
FFLAGS += $(INC)

$(EXEC): Driver.f90 $(MODS) $(SUBS)
	$(F90) $(FFLAGS) Driver.f90 $(MODS) $(SUBS) -o $(@) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(MODS) *.mod

