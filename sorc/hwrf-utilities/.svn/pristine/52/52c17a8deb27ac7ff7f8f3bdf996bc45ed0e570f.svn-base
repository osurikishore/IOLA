
include ../../configure.hwrf

LIBS = -L../../libs -lg2 -lw3nco_i4r4 -lbacio_4 $(LDFLAGS_GRIB2_LIB)
EXE=../../exec/satgrib2.exe

all: $(EXE)

clean:
	rm -f *.o *.mod satgrib2 a.out $(EXE) *~

$(EXE): satgrib2.F90 gds2gdt.f
	$(SFC) $(FFLAGS)      \
		 -o $(EXE) satgrib2.F90 gds2gdt.f $(LIBS)
