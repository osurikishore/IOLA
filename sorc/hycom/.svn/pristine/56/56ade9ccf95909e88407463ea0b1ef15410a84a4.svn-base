include ../../configure.hycom
FC=$(SFC)

SHELL         = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .F .f .o


LIBS+=   -L../../../hwrf-utilities/libs/ -lw3nco_i4r4 -lbacio
#LIBS+=   $(W3NCO_LIB4) $(BACIO_LIB4)
CPPFLAGS =  -P -D$(BYTE_ORDER) -DWRF -DLINUX -DSERIAL_IO -I$(WRF_DIR)/inc



MODS = mod_xc.o mod_za.o wtime.o
EXEC     =      ../../exec

program=  ../../exec/ofs_latlon


CMD=ofs_latlon

SRCS=	grib.f

OBJS=	${SRCS:.f=.o}


all:	$(program)

clean:
	/bin/rm -f *.o *.a *.mod *.F.f $(program)

$(program):	$(MODS) $(OBJS) zh.o 
	${FC} ${FCFLAGS} -o $(@) $(OBJS) $(MODS) zh.o $(LIBS) 

grib.o: grib.f
mod_xc.o:   mod_xc.F
mod_za.o:   mod_za.F mod_xc.o
zh.o:       zh.F 


