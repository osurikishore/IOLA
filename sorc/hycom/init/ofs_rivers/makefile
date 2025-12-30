include ../../configure.hycom

FC=$(SFC)

#LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r4 -lbacio
LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

INC = ../../libs/src/ofs_mods
FFLAGS+= -I$(INC)

MODS=      mod_xc.o mod_za.o zh.o wtime.o mod_zb.o
OBJ  =     rivers.o 
ALLOBJ = $(MODS) $(OBJ)

CMD	= ../../exec/ofs_rivers

$(CMD): $(ALLOBJ)
	$(FC) $(LDFLAGS) -o $@ $(ALLOBJ) $(LIBS)

.PHONY: clean

clean:
	$(RM) $(ALLOBJ) $(CMD) *.mod *.o

