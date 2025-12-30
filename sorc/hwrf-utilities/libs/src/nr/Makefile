
LIB = ../../libnr.a

include ../../../configure.hwrf

MODS = 

SUBS = *.o 

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(LIB)

.SUFFIXES: .f .o

.f.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.f  
