#--------------------------------------------------------------------
# This makefile is for the ip library
#
# The makefile expects two variables set from the command line
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB = $(OUTPUT_LIB)

#
# configuration based on architecture and compiler
include ../../../configure.hwrf

FC=$(SFC)
CC=$(SCC)
FFLAGS=$(FFLAGS_BUFR) $(FFLAGS_MORE) -I../../mods/sigio -DLINUX $(OMP_FFLAGS)
#FFLAGS+=-O0 -g -traceback -check all
CFLAGS=$(CFLAGS_BUFR) $(CFLAGS_MORE) -DLINUX

#
# Object files for this library - all are f90
OBJS_F   =
OBJS_F77 =                                    \
  gausslat.o gcdist.o gdsawt.o gdswiz00.o     \
  gdswiz01.o gdswiz03.o gdswiz04.o gdswiz05.o \
  gdswizc9.o gdswizca.o gdswizcb.o gdswizcd.o \
  gdswiz.o gdswzd00.o gdswzd01.o gdswzd03.o   \
  gdswzd04.o gdswzd05.o gdswzdc9.o gdswzdca.o \
  gdswzdcb.o gdswzdcd.o gdswzd.o ijkgds0.o    \
  ijkgds1.o ijkgds.o ipmerge2.o ipolates.o    \
  ipolatev.o ipsector.o ipspaste.o ipxetas.o  \
  ipxwafs2.o ipxwafs3.o ipxwafs.o makgds.o    \
  movect.o polateg0.o polateg1.o polateg4.o   \
  polates0.o polates1.o polates2.o polates3.o \
  polates4.o polates6.o polatev0.o polatev1.o \
  polatev2.o polatev3.o polatev4.o polatev6.o \
  polfixs.o polfixv.o
OBJS_C   =
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

#
# Source files
SRCS_F   = $(OBJS_F:.o=.f90)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETS
all:
	OUTPUT_LIB=../../libip_i4r8.a \
	FFLAGS_MORE="$(FFLAGS_DOUBLE)" \
	CFLAGS_MORE="" \
	make -f Makefile clean lib
	FFLAGS_MORE="$(FFLAGS_SINGLE)" \
	CFLAGS_MORE="" \
	OUTPUT_LIB=../../libip_i4r4.a \
	make -f Makefile clean lib

#---------------------
# Main compilation and archiving of the library

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $?

clean:
	$(RM) $(OBJS)

bare: clean
	$(RM) ../../libip*a

# Override the .F to be fixed-form pre-processed
.F.o:
	$(RM) $*.mod $*.o
	$(FC) $(FPPFLAGS) $(FFLAGS) -c $*.F

.f.o:
	$(RM) $*.mod $*.o
	$(FC) $(FFLAGS) -c $*.f 

.IGNORE:
.PHONY: clean
