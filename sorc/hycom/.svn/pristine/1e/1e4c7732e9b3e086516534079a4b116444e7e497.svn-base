include ../../configure.hycom
#F90=$(SF90)
#FC=$(SF90)

LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r4 -lbacio -lg2 -L/usr/lib64 -ljasper
INC += ../../../hwrf-utilities/libs/mods/g2
FFLAGS+= -I$(INC)

OBJS    = getgds2.o
CMD    = ../../exec/ofs_getgds2

$(CMD):        $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(CMD)

