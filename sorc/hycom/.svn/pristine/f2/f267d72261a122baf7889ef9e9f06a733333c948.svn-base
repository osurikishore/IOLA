include ../../configure.hycom
F90=$(SF90)

LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3nco_i4r4 -lbacio
#LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

OBJS    = getkpds.o
CMD    = ../../exec/ofs_getkpds

$(CMD):        $(OBJS)
	$(SFC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(CMD)

