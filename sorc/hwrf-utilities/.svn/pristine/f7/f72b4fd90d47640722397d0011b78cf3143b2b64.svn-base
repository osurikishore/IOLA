include ../../../configure.hwrf

OBJS = bacio.o byteswap.o baciof.o bafrio.o chk_endianc.o

FFLAGS = $(FFLAGS_BUFR) $(FFLAGS_MORE) -DLINUX
CFLAGS+= -W

all:
	LIB1=../../libbacio_8.a \
	FFLAGS_MORE="$(FFLAGS_DOUBLE2)" \
	CFLAGS_MORE="" \
	INCMOD=../../mods/bacio_8 \
	make -f Makefile clean lib1
	FFLAGS_MORE="$(FFLAGS_SINGLE)" \
	CFLAGS_MORE="" \
	LIB1=../../libbacio.a \
	LIB2=../../libbacio_4.a \
	INCMOD=../../mods/bacio_4 \
	make -f Makefile clean lib1 lib2

lib1: $(LIB1)
lib2: $(LIB2)

$(LIB1): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB1) $(OBJS)
	$(MKDIR) $(INCMOD)
	$(CP) *.mod $(INCMOD)

$(LIB2): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB2) $(OBJS)
	$(MKDIR) $(INCMOD)
	$(CP) *.mod $(INCMOD)

.f.o:
	$(SFC) $(FFLAGS) -c $<
.c.o:
	$(SCC) $(CFLAGS) $(DEFS) -c $<
.PHONY: clean

clean:
	$(RM) *.o *.mod *.a

