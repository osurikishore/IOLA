include ../../configure.hycom

FC=$(SFC)

LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3nco_i4r4 -lbacio -lbufr_i4r4
#LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4) $(BUFR_LIB4)

OBJS    = getdata_river.o
CMD    = ../../exec/ofs_getdata_river

$(CMD):        $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(CMD)

