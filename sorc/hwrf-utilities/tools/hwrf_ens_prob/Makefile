include ../../configure.hwrf

PROBOBJS = hwrf_ens_prob.o ens_prob_module.o ens_output_module.o
PROBEXEC=../../exec/hwrf_ens_prob.exe

ATCFOBJS = atcf_prob_module.o atcf_prob.o ens_prob_module.o ens_output_module.o
ATCFEXEC=../../exec/hwrf_atcf_prob.exe

EXECS=$(PROBEXEC) $(ATCFEXEC)

LIBDIR=../../libs
LIBSRC=../../libs/mods/hwrfutil
LIBS = -L$(LIBDIR) -lhwrfutil_i4r4 -lnetcdf
LIB = $(LIBDIR)/libhwrfutil_i4r4.a
FFLAGS += $(NETCDF_FFLAGS) -I$(LIBSRC) $(OMP_FFLAGS)
LDFLAGS += $(NETCDF_LDFLAGS) $(OMP_FFLAGS)
F90=$(SFC) 
CC=$(SCC)

.SUFFIXES: # GNU Make bug workaround
.SUFFIXES: .f90 .o .F90

# Lines from here on down should not need to be changed.
all: $(EXECS)

bare: clean
	$(RM) $(EXECS) *.a

clean:
	$(RM) $(OBJS) *.mod *~ a.out *.o

$(PROBEXEC): $(PROBOBJS)  $(LIB)
	$(F90) -o $(@) $(PROBOBJS)  $(LIBS) $(LDFLAGS) $(FFLAGS)
$(ATCFEXEC): $(ATCFOBJS)  $(LIB)
	$(F90) -o $(@) $(ATCFOBJS)  $(LIBS) $(LDFLAGS) $(FFLAGS)

.PHONY: clean bare

.F90.o:
	$(F90) $(FFLAGS) -c -o $*.o $*.F90

atcf_prob_module.o: ens_prob_module.o ens_output_module.o $(LIB)
atcf_prob.o: atcf_prob_module.o $(LIB)
glob_test.o: $(LIB)

ens_prob_module.o: $(LIB)
ens_output_module.o: ens_prob_module.o $(LIB)
hwrf_ens_prob.o: ens_prob_module.o ens_output_module.o $(LIB)

