include ../../configure.hwrf

SRCS=	gds2gdt.f pds2pdt.f pds2pdtens.f cnvgrib.f cnv12.f cnv22.f \
        cnv21.f gdt2gds.f makepds.f putgbexn.f makepdsens.f setbit.f
OBJS=	gds2gdt.o pds2pdt.o pds2pdtens.o cnvgrib.o cnv12.o cnv22.o \
        cnv21.o gdt2gds.o makepds.o putgbexn.o makepdsens.o setbit.o
FC =		$(SFC)
FFLAGS = $(FFLAGS_BUFR) -I../../libs/mods/g2
LIBS = -L../../libs -lg2 -lw3nco_i4r4 -lbacio_4 $(LDFLAGS_GRIB2_LIB)
CMD =		../../exec/cnvgrib.exe

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
