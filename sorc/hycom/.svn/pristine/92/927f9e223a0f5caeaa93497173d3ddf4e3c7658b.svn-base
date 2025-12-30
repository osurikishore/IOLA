include ../../configure.hycom
#
# If it is invoked by the command line
#	make -f makefile
#
FC            = $(SFC)
CMD =      ../../exec/hwrf_select_domain
EXEC =     ../../exec

OBJS=        select_domain.o

OBJST=

LIBS =        

$(CMD):      $(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS)

#$(OBJS):    
#	$(SFC) $(FFLAGS) -c $<

#$(OBJST):    
#	$(FC) $(FFLAGST) -c $<

clean :
	-rm -f *.o $(CMD)

install :
	-mv $(CMD) $(EXEC)


