
LIB = ../../libkdtree.a
INCMOD = ../../mods/kdtree

include ../../../configure.hwrf

MODS = kdtree2.o

SUBS = 

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(MKDIR) $(INCMOD)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	mv *.mod $(INCMOD)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(INCMOD)/*.mod $(LIB)

.SUFFIXES: .F90 .o

.F90.o:
	$(SFC) $(FFLAGS) -c $*.F90
