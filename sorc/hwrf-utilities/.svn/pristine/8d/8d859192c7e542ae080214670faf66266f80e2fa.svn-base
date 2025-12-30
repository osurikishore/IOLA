
LIB = ../../libslint.a
INCMOD = ../../mods/slint

include ../../../configure.hwrf

MODS = kd.o uvxyz.o slint.o

SUBS = 

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(MKDIR) $(INCMOD)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	mv *.mod $(INCMOD)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(INCMOD)/*.mod $(LIB)

.SUFFIXES: .F90 .f .o

.F90.o:
	$(SFC) $(FFLAGS) -c $*.F90

.f.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.f  
