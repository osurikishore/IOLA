include ../../configure.hycom
FC=$(SFC)
#
# If it is invoked by the command line
#	make -f makefile
#
CMD =      ../../exec/hwrf_get_layers
EXEC = 	../../exec

#FC = ncepxlf90
#FX = ncepxlf

OBJS=        get_layers.o

OBJST=       

LIBS =        

all: 	$(CMD)

$(CMD):      $(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS)

install:
	-mv $(CMD) $(EXEC)

clean :
	-rm -f *.o $(CMD)
      
