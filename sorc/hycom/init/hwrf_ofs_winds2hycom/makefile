include ../../configure.hycom
F90=$(SF90)
FC=$(SF90)

SRCS=  mod_gribio.f90  mod_flags.f90 constants.f90 horiz_interp.f90 mod_geom.f90  wind2hycom.f90  flush.f90

OBJS=	${SRCS:.f90=.o}

EXEC= ../../exec

#LIBS    += -L../../libs -lofs_mods -L$(LIB_W3_PATH) -lw3_i4r8 -L$(LIB_BACIO_PATH) -lbacio
#LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r8 -lbacio
LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3nco_i4r4 -lbacio
#LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

INC = ../../libs/src/ofs_mods
FFLAGS+= -I$(INC)
#CPPFLAGS+= -DIFC

CMD   =		../../exec/hwrf_ofs_wind2hycom

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build swathps.

.IGNORE:
.SUFFIXES: .f .o .F .F90 .f90

all:	$(CMD)

install:
	-mv $(CMD) $(EXEC)

$(CMD):		$(OBJS) $(OBJS_f) 
	$(F90) $(LDFLAGS) -o $(@) $(OBJS) $(OBJS_f) $(LIBS)

clean:
	-rm -f $(OBJS) $(OBJS_f) *.mod $(CMD)

allclean:	clean all

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile

test_ud: test_ud.o
	$(SFC) $(LDFLAGS) -o test_ud test_ud.o

test_ud.o:	
	$(SFC) $(FFLAGS) -c test_ud.f90 -o test_ud.o
