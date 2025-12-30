include ../../configure.hycom
FC=$(SFC)
F90=$(SF90)
 
SRCS=   correct_wind.f90 
SRCS_f= get_abfld1.f       

OBJS=	${SRCS:.f90=.o}
OBJS_f=	${SRCS_f:.f=.o}

EXEC=../../exec

#LIBS    += -L../../libs -lofs_mods -L$(LIB_W3_PATH) -lw3_i4r8 -L$(LIB_BACIO_PATH) -lbacio
#LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r8 -lbacio
#LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r4 -lbacio
LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

INC = ../../libs/src/ofs_mods
FFLAGS+= -I$(INC)

CMD   =		../../exec/hwrf_ofs_correct_wind

.IGNORE:
.SUFFIXES: .c .f .F .F90 .f90 .o
all:	$(CMD)

$(CMD):		$(OBJS) $(OBJS_f) 
	$(SF90) $(LDFLAGS) -o $(@) $(OBJS) $(OBJS_f) $(LIBS)

clean:
	$(RM) -f *.o *.mod $(CMD)
