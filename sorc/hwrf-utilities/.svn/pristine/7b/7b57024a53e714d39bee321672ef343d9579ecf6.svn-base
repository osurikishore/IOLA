#!/bin/csh

if ( -d /lfs4 ) then
    setenv platform  "jet"
else if ( -d  /scratch1 ) then
    setenv platform  "hera"
else if ( -d /gpfs/hps ) then
    setenv platform  "wcosscray"
else if ( -d /work/noaa ) then
    setenv platform  "orion"
else
    setenv platform  "unknown"
endif

if ( ${platform} == "jet" ) then
     setenv NCEPLIBS /mnt/lfs4/HFIP/hfv3gfs/gwv/l0716/lib/
     setenv BUFR_LIB4 ${NCEPLIBS}/bufr/v11.3.1/intel/libbufr_v11.3.1_4_64.a
     setenv BUFR_LIB8 ${NCEPLIBS}/bufr/v11.3.1/intel/libbufr_v11.3.1_8_64.a
     setenv BUFR_LIBd ${NCEPLIBS}/bufr/v11.3.1/intel/libbufr_v11.3.1_d_64.a
     module use /mnt/lfs4/HFIP/hfv3gfs/gwv/l0716/lib/modulefiles
     module load  bufr-intel-sandybridge/11.3.1
     /bin/cp ${BUFR_LIB4} ../../libbufr_i4r4.a
     /bin/cp ${BUFR_LIB8} ../../libbufr_i8r8.a 
     /bin/cp ${BUFR_LIBd} ../../libbufr_i4r8.a
else if ( ${platform} == "hera" ) then
     setenv NCEPLIBS /scratch2/NCEPDEV/nwprod/NCEPLIBS/
     setenv BUFR_LIB4 ${NCEPLIBS}/compilers/intel/18.0.5.274/lib/libbufr_v11.3.0_4_64.a
     setenv BUFR_LIB8 ${NCEPLIBS}/compilers/intel/18.0.5.274/lib/libbufr_v11.3.0_8_64.a
     setenv BUFR_LIBd ${NCEPLIBS}/compilers/intel/18.0.5.274/lib/libbufr_v11.3.0_d_64.a
     module use /scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles 
     module load bufr/11.3.0
     /bin/cp ${BUFR_LIB4} ../../libbufr_i4r4.a
     /bin/cp ${BUFR_LIB8} ../../libbufr_i8r8.a 
     /bin/cp ${BUFR_LIBd} ../../libbufr_i4r8.a
else if ( ${platform} == "wcosscray" ) then
     echo "Need to Add Here"
else if ( ${platform} == "orion" ) then
     setenv NCEPLIBS  /apps/contrib/NCEPLIBS/orion/lib/
     setenv BUFR_LIB4 ${NCEPLIBS}/libbufr_v11.3.0_4_64.a
     setenv BUFR_LIB8 ${NCEPLIBS}/libbufr_v11.3.0_8_64.a
     setenv BUFR_LIBd ${NCEPLIBS}/libbufr_v11.3.0_d_64.a
     module use /apps/contrib/NCEPLIBS/orion/modulefiles
     module load  bufr/11.3.0
     /bin/cp ${BUFR_LIB4} ../../libbufr_i4r4.a
     /bin/cp ${BUFR_LIB8} ../../libbufr_i8r8.a
     /bin/cp ${BUFR_LIBd} ../../libbufr_i4r8.a
else
     echo "Unknown Platform  [Not Jet, Hera, Wcosscray or Orion (~ _~;)]"
endif

