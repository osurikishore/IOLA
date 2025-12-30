include ../../configure.hycom
#
# If it is invoked by the command line
#	make -f makefile
#
CMD =      ../../exec/hwrf_ofs_fwind
EXEC = 	../../exec

FC = $(SFC)
FX = $(SFC)

OBJS=        fwind.o

OBJST=       

LIBS =        

all: 	$(CMD)

$(CMD):      $(OBJST) $(OBJS)
	$(FX) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS)

#$(OBJS):    
#	$(FX) $(FFLAGS) -c $<
#
#$(OBJST):    
#	$(FX) $(FFLAGST) -c $<

install:
	-mv $(CMD) $(EXEC)

clean :
	-rm -f *.o $(CMD)
      
