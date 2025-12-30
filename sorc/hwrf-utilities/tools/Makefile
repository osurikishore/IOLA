# Must include configure.hwrf here or CPPFLAGS and CPP are overrided
# as "" on AIX.
include ../configure.hwrf

SHELL=/bin/sh

# HWRF tools makefile

SUBDIRS = hwrf_tempdrop hwrf_obstobufr hwrf_regrid_merge wgrib2 \
        satgrib2 cnvgrib copygb hwrf_wrfbdy_tinterp hwrf_ens_prob \
	grbindex grb2index hwrf_atcf_to_stats hwrf_aux_rw hwrf_prepbufr  \
	hwrf_gridgenfine hwrf_prep_hybrid hwrf_read_indi  \
	hwrf_supvit hwrf_swath wgrib ndate mdate nhour hwrf_bin_io \
	hwrf_afos hwrf_combinetrack serpoe hwrf_wrfout_newtime \
	hwrf_htcfstats hwrf_readtdrinfo \
	ships hwrf_bdy_update hwrf_netcdf_grads hwrf_binary_grads \
	bufr_remorest hwrf_nhc_products mpiserial \
	hwrf_blend_gsi mpi_example hwrf_metgrid_levels \
	hwrf_final_merge hwrf_ensemble hwrf_da_update_analysis

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; make ); \
	done

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
	   ( cd $$dir; echo "Making $@ in `pwd`" ; \
	     make $@) ; \
	done
