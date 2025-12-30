SHELL=/bin/sh

include ../../../configure.hwrf

OBJS_MODA = \
moda_bitbuf.o moda_bufrmg.o moda_bufrsr.o  \
moda_comprs.o moda_comprx.o moda_h4wlc.o moda_idrdm.o moda_ifopbf.o  \
moda_ival.o moda_ivttmp.o moda_lushr.o moda_mgwa.o moda_mgwb.o  \
moda_msgcwd.o moda_msglim.o moda_msgmem.o moda_mstabs.o moda_nmikrp.o  \
moda_nulbfr.o moda_rdmtb.o moda_s01cm.o moda_sc3bfr.o moda_stbfr.o  \
moda_stcode.o moda_tababd.o moda_tables.o moda_ufbcpl.o moda_unptyp.o  \
moda_usrbit.o moda_usrint.o moda_usrtmp.o moda_xtab.o

OBJS_MODV= modv_MAXCD.o  \
modv_MAXJL.o modv_MAXMEM.o modv_MAXMSG.o modv_MAXSS.o modv_MAXTBA.o  \
modv_MAXTBB.o modv_MAXTBD.o modv_MXCDV.o modv_MXCSB.o modv_MXDXTS.o  \
modv_MXH4WLC.o modv_MXLCC.o modv_MXMSGL.o modv_MXMTBB.o modv_MXMTBD.o  \
modv_MXS01V.o modv_NFILES.o

OBJS_NOTMOD = \
arallocc.o ardllocc.o bort_exit.o ccbfl.o cmpia.o cobfl.o cpmstabs.o  \
crbmg.o cread.o cwbmg.o icvidx.o nummtb.o rbytes.o restd.o stseq.o  \
wrdesc.o adn30.o atrcpt.o bfrini.o blocks.o bort2.o bort.o bvers.o  \
cadn30.o capit.o chekstab.o chrtrna.o chrtrn.o cktaba.o closbf.o  \
closmg.o cmpmsg.o cmsgini.o cnved4.o conwin.o copybf.o copymg.o  \
copysb.o cpbfdx.o cpdxmm.o cpymem.o cpyupd.o datebf.o datelen.o digit.o  \
drfini.o drstpl.o dumpbf.o dxdump.o dxinit.o dxmini.o elemdx.o errwrt.o  \
fstag.o getabdb.o getbmiss.o getlens.o getntbe.o gets1loc.o gettagpr.o  \
gettbh.o getvalnb.o getwin.o hold4wlc.o i4dy.o ibfms.o icbfms.o  \
ichkstr.o icmpdx.o icopysb.o idn30.o idxmsg.o ifbget.o ifxy.o  \
igetdate.o igetfxy.o igetntbi.o igetntbl.o igetprm.o igetsc.o igettdi.o  \
inctab.o invcon.o invmrg.o invtag.o invwin.o iok2cpy.o iokoper.o ipkm.o  \
ipks.o ireadmg.o ireadmm.o ireadns.o ireadsb.o ishrdx.o isize.o  \
istdesc.o iupb.o iupbs01.o iupbs3.o iupm.o iupvs01.o jstchr.o jstnum.o  \
lcmgdf.o lmsg.o lstjpb.o makestab.o maxout.o mesgbc.o mesgbf.o minimg.o  \
mrginv.o msgfull.o msgini.o msgupd.o msgwrt.o mtinfo.o mvb.o nemdefs.o  \
nemock.o nemspecs.o nemtab.o nemtba.o nemtbax.o nemtbb.o nemtbd.o  \
nenubd.o nevn.o newwin.o nmsub.o nmwrd.o numbck.o numtab.o numtbd.o  \
nvnwin.o nwords.o nxtwin.o openbt.o openmb.o openmg.o pad.o padmsg.o  \
parstr.o parusr.o parutg.o pkb.o pkbs1.o pkc.o pkftbv.o pktdd.o pkx.o  \
posapx.o rcstpl.o rdbfdx.o rdcmps.o rdmemm.o rdmems.o rdmgsb.o rdmsgw.o  \
rdmtbb.o rdmtbd.o rdtree.o rdusdx.o readdx.o readerme.o readlc.o  \
readmg.o readmm.o readns.o reads3.o readsb.o rewnbf.o rjust.o rsvfvm.o  \
rtrcptb.o rtrcpt.o seqsdx.o setblock.o setbmiss.o sntbbe.o sntbde.o  \
status.o stbfdx.o stdmsg.o stndrd.o stntbia.o stntbi.o strcln.o  \
strcpt.o string.o strnum.o strsuc.o tabent.o tabsub.o trybump.o  \
ufbcnt.o ufbcpy.o ufbcup.o ufbdmp.o ufbevn.o ufbget.o ufbin3.o ufbint.o  \
ufbinx.o ufbmem.o ufbmex.o ufbmms.o ufbmns.o ufbovr.o ufbpos.o ufbqcd.o  \
ufbqcp.o ufbrep.o ufbrms.o ufbrp.o ufbrw.o ufbseq.o ufbsp.o ufbstp.o  \
ufbtab.o ufbtam.o ufdump.o upbb.o upb.o upc.o upds3.o upftbv.o ups.o  \
uptdd.o usrtpl.o valx.o wrcmps.o wrdxtb.o writcp.o writdx.o writlc.o  \
writsa.o writsb.o wrtree.o wtstat.o arallocf.o ardllocf.o exitbufr.o  \
irev.o isetprm.o openbf.o pkvs01.o readmt.o wrdlen.o 

OBJS=$(OBJS_MODA) $(OBJS_MODV) $(OBJS_NOTMOD)

PRM = bufrlib.prm

FC=$(SFC)
CC=$(SCC)
CFLAGS=$(CFLAGS_BUFR) $(CFLAGS_BUFR_MORE) -DUNDERSCORE -D$(BYTE_ORDER) -DDYNAMIC_ALLOCATION -DLINUX
FPPFLAGS =  $(FPPFLAGS_BUFR_MORE) -DDYNAMIC_ALLOCATION -DLINUX
FFLAGS=$(FFLAGS_BUFR) $(FFLAGS_BUFR_MORE) $(FPPFLAGS) -I../../mods/sigio/ -DLINUX
LIB=$(OUTPUT_LIB)

.F.o:
	$(SFC) $(FFLAGS) -c $< -o $@

.f.o:
	$(SFC) $(FFLAGS) -c $< -o $@

.c.o:
	$(SCC) $(CFLAGS) -c $< -o $@

all: i4r8 i4r4

i4r8: 
	set -xue ; \
	export OUTPUT_LIB=../../libbufr_i4r8.a ; \
	export OUTPUT_MOD=../../mods/bufr_i4r8 ; \
	export FFLAGS_BUFR_MORE="$(FFLAGS_DOUBLE)" ; \
	export FPPFLAGS_BUFR_MORE="-DSUPERSIZE_BUILD" ; \
	export CFLAGS_BUFR_MORE="-DSUPERSIZE_BUILD" ; \
	make -f Makefile clean bufrlib.prm ; \
	cflags_defs='' ; \
        for bprm in MAXNC MXNAF ; do \
            bprmval=`grep " $${bprm} = " bufrlib.prm | cut -f2 -d= | cut -f2 -d" "` ; \
            cflags_defs="$${cflags_defs} -D$${bprm}=$${bprmval}" ; \
        done ; \
        for gvar in NFILES MAXCD ; do \
            gvarval=`$(FPP) $(FPPFLAGS) $$FPPFLAGS_BUFR_MORE modv_$$gvar.F | grep " $${gvar} = " | cut -f2 -d= | cut -f2 -d" "` ; \
            cflags_defs="$${cflags_defs} -D$${gvar}=$${gvarval}" ; \
        done ; \
	echo "cflags_defs = $${cflags_defs}" ; \
	OUTPUT_LIB=../../libbufr_i4r8.a \
	FFLAGS_BUFR_MORE="$(FFLAGS_DOUBLE)" \
	FPPFLAGS_BUFR_MORE="-DSUPERSIZE_BUILD" \
	CFLAGS_BUFR_MORE="-DSUPERSIZE_BUILD $$cflags_defs" \
	make -f Makefile lib

i4r4: 
	set -xue ; \
	export OUTPUT_LIB=../../libbufr_i4r4.a ; \
	export OUTPUT_MOD=../../mods/bufr_i4r8 ; \
	export FFLAGS_BUFR_MORE="$(FFLAGS_SINGLE)" ; \
	export FPPFLAGS_BUFR_MORE="-DSUPERSIZE_BUILD" ; \
	export CFLAGS_BUFR_MORE="-DNORMAL_BUILD" ; \
	make -f Makefile clean bufrlib.prm ; \
	cflags_defs='' ; \
        for bprm in MAXNC MXNAF ; do \
            bprmval=`grep " $${bprm} = " bufrlib.prm | cut -f2 -d= | cut -f2 -d" "` ; \
            cflags_defs="$${cflags_defs} -D$${bprm}=$${bprmval}" ; \
        done ; \
        for gvar in NFILES MAXCD ; do \
            gvarval=`$(FPP) $(FPPFLAGS) $$FPPFLAGS_BUFR_MORE modv_$$gvar.F | grep " $${gvar} = " | cut -f2 -d= | cut -f2 -d" "` ; \
            cflags_defs="$${cflags_defs} -D$${gvar}=$${gvarval}" ; \
        done ; \
	echo "cflags_defs = $${cflags_defs}" ; \
	OUTPUT_LIB=../../libbufr_i4r4.a \
	FFLAGS_BUFR_MORE="$(FFLAGS_DOUBLE)" \
	FPPFLAGS_BUFR_MORE="-DNORMAL_BUILD" \
	CFLAGS_BUFR_MORE="-DNORMAL_BUILD $$cflags_defs" \
	make -f Makefile lib

$(OBJS_MODA): $(OBJS_MODV)
$(OBJS_NOTMOD): $(OBJS_MODA) $(OBJS_MODV)

lib:	$(PRM) $(LIB)

$(LIB):	$(OBJS)
	$(AR) -ruv $(LIB) $(OBJS)
	$(MKDIR) $(OUTPUT_MOD)
	$(CP) *.mod $(OUTPUT_MOD)

bufrlib.prm: bufrlib.PRM
	$(FPP) $(FPPFLAGS) bufrlib.PRM > bufrlib.prm

clean:
	rm -f *.o *.mod $(LIB) bufrlib.prm \
	moda_idrdm.f moda_lushr.f moda_ival.f moda_ivttmp.f


# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

$(OBJS):	bufrlib.prm
bufrlib.prm:	bufrlib.PRM
