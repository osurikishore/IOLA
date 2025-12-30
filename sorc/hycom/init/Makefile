# Must include configure.hycom here or CPPFLAGS and CPP are overrided
# as "" on AIX.
include ../configure.hycom

SHELL=/bin/sh
#HYCOM init makefile

SUBDIRS = hwrf_get_rtofs \
	hwrf_archv2restart \
	hwrf_gfs2ofs2 \
	hwrf_ofs_correct_wind \
	hwrf_ofs_fwind \
	hwrf_ofs_winds2hycom \
	hwrf_restart2restart \
	hwrf_select_domain \
	hwrf_subregion \
	hwrf_ofs_timeinterp_forcing \
	hwrf_rtofs_edit_gfsfrc \
	ofs_correct_forcing \
	ofs_getdata_river  \
	ofs_getkpds  \
	ofs_qc_river  \
	ofs_rivers   \
	ofs_forcing_add \
	hwrf_hycom2raw \
	iola_ofs_timeinterp_forcing

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; make ); \
	done

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; \
	   make $@) ; \
	done

