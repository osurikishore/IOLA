################################################################################
include ../../configure.hwrf
SHELL=	/bin/sh
CMD=	../../exec/hwrf_regrid_merge.exe
LOPTS=  $(OMP_FFLAGS)
INCS=	
OBJS=	module_regrid_merge.o hwrf_regrid_merge.o from_copygb.o
LIBS=  -L../../libs -lw3nco_i4r8 -lbacio_4  -lip_i4r8 -lsp_i4r8 -lw3nco_i4r8

FC=$(SFC)
FFLAGS=$(FFLAGS_BUFR) $(FFLAGS_DOUBLE) $(OMP_FFLAGS)
#FFLAGS+=-warn -g -traceback -O0 -check all
################################################################################
all: $(CMD)
$(CMD):	$(OBJS)
	$(SFC) $(LOPTS) $(OBJS) $(LIBS) -o $(CMD)
clean:	
	rm -f $(OBJS)
.SUFFIXES: # Workaround for GNU Make bug
.SUFFIXES: .F90 .o .exe .F
.PHONY: clean
.F90.o:
	$(RM) $@ $*.mod
	$(F90) $(FFLAGS) -c $*.F90 
.F.o:
	$(RM) $@ $*.mod
	$(F90) $(FFLAGS) -c $*.F
