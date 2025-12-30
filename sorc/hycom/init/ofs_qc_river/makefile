include ../../configure.hycom

FC=$(SFC)

#LIBS    += -L../../../hwrf-utilities/libs/ -lw3_i4r4 -lbacio
LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

OBJS    = qc_river.o
CMD    = ../../exec/ofs_qc_river

$(CMD):        $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) $(CMD)

