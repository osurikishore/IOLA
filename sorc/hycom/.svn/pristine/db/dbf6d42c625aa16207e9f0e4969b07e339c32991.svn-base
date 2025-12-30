# Must include configure.hycom here or CPPFLAGS and CPP are overrided
# as "" on AIX.
include ../configure.hycom

SHELL=/bin/sh
#HYCOM post makefile

SUBDIRS = hwrf_ab2data \
	hwrf_get_layers \
	ofs_latlon

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; make ); \
	done

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; \
	   make $@) ; \
	done

