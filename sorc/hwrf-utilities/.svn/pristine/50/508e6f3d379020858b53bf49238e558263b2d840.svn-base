
LIB = ../../libfson.a
INCMOD = ../../mods/fson

include ../../../configure.hwrf

MODS = fson_string_m.o fson_value_m.o fson_path_m.o fson.o

SUBS = 

OBJS = $(MODS) $(SUBS)

$(LIB): $(OBJS)
	$(MKDIR) $(INCMOD)
	$(AR) $(ARFLAGS) $(AFLAGS) $@ $(OBJS)
	mv *.mod $(INCMOD)
	$(RM) -f *.o

clean:
	$(RM) -f *.o $(INCMOD)/*.mod $(LIB)

.SUFFIXES: .f90 .f .o

.f90.o:
	$(SFC) $(FFLAGS) -c $*.f90

.f.o:
	$(SFC) $(FFLAGS) $(FIXED) -c $*.f  
