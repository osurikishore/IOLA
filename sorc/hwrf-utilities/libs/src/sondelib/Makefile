
LIB = ../../libsondelib.a
INCMOD = ../../mods/sondelib

include ../../../configure.hwrf

MODS = aoml_routines.o

SUBS = thermo_subs.o \
       sondelib.o

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(MKDIR) $(INCMOD)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(INCMOD)/*.mod $(LIB)

.SUFFIXES: .F90 .f .o

.F90.o:
	$(SFC) $(FFLAGS) -c $*.F90

.f.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.f  
