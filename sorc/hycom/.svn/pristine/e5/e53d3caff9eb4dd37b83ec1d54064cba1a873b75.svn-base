include ../../configure.hycom
F90=$(SF90)
FC=$(SFC)

#LIBS    += -L../../libs -lofs_mods -L../../../hwrf-utilities/libs/ -lw3_i4r4 -lbacio 
LIBS    += -L../../libs -lofs_mods $(W3NCO_LIB4) $(BACIO_LIB4)

INC = ../../libs/src/ofs_mods
FFLAGS+= -I$(INC)

MODS =      mod_xc.o mod_za.o zh.o wtime.o
OBJ  =     get_abfld1.o
ALLOBJ = $(MODS) $(OBJ)

.SUFFIXES: .c .f .F .F90 .f90 .o

CMD	= ../../exec/ofs_correct_forcing

$(CMD): $(ALLOBJ)  correct_forcing.o
	$(DM_FC) $(LDFLAGS) correct_forcing.o $(ALLOBJ) $(LIBS) -o $@

.PHONY: clean

clean:
	$(RM) $(ALLOBJ) $(CMD) *.mod *.o


