
LIB = ../../libspline.a
INCMOD = ../../mods/spline

include ../../../configure.hwrf

MODS = spline.o

SUBS = 

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(MKDIR) $(INCMOD)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(INCMOD)/*.mod $(LIB)

.SUFFIXES: .f90 .f .o

.f90.o:
	$(SFC) $(FFLAGS) -c $*.f90

.f.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.f  
