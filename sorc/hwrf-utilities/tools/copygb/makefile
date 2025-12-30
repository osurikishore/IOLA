################################################################################
include ../../configure.hwrf
SHELL=	/bin/sh
CMD=	../../exec/copygb.exe
FOPTS=	$(FFLAGS) $(FFLAGS_DOUBLE) $(OMP_FFLAGS)
LOPTS=  $(OMP_FFLAGS)
INCS=	
OBJS=	copygb.o
LIBS=  -L../../libs -lw3nco_i4r8 -lbacio_4  -lip_i4r8 -lsp_i4r8
################################################################################
$(CMD):	$(OBJS)
	$(SFC) $(LOPTS) $(OBJS) $(LIBS) -o $(CMD)
$(OBJS):	$(INCS)
copygb.o:	copygb.F
	$(SFC) -c $(FOPTS) $*.F
clean:	
	-rm -f $(OBJS)
.PHONY: clean
