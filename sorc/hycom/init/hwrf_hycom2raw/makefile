include ../../configure.hycom

FC=$(SFC)
#
# If it is invoked by the command line
#	make -f makefile
#
CMD =      ../../exec/hwrf_hycom2raw
EXEC = 	../../exec


OBJS=        hycom2raw.o hycom_endian_io.o parse.o

OBJST=       

LIBS =        

.SUFFIXES: .c .F .f .o
#
# rules.
#

all: 	$(CMD)

$(CMD):      $(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS)

install:
	-mv $(CMD) $(EXEC)

clean :
	-rm -f *.o $(CMD)
      
