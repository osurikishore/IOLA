include ../../configure.hwrf


OBJS =	grb2index.o

LIBS  = -L../../libs -lg2 -lw3nco_i4r4 -lbacio_4 -lw3_i4r4 
EXEC  =  ../../exec/grb2index.exe

# Lines from here on down should not need to be changed.

$(EXEC): $(OBJS)
	$(SFC) $(LDFLAGS) $(FFLAGS) -o $(@) grb2index.f $(LIBS)

.PHONY: clean

clean:
	$(RM) $(OBJS) *.mod 
