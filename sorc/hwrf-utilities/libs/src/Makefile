SHELL=/bin/sh
# HWRF LIBS makefile

.PHONY: all clean

all: make_fson make_slint make_sondelib make_slatec make_nr \
     make_bacio make_blas make_bufr make_sfcio make_sp make_w3nco make_w3emc \
     make_sigio make_hwrfutil make_g2 make_g2tmpl make_ip make_g2c make_gctpc \
     make_nemsio make_kdtree make_spline

clean: clean_fson clean_slint clean_sondelib clean_slatec clean_nr \
       clean_bacio clean_blas clean_bufr clean_sfcio clean_sp clean_w3nco \
       clean_w3emc clean_sigio clean_hwrfutil clean_g2 clean_g2tmpl \
       clean_ip clean_g2c clean_gctpc clean_nemsio clean_kdtree clean_spline

bare: clean
	rm -f ../*.a
	rm -rf ../mods

make_spline:
	rm -f make_spline
	set -e; cd spline; make -j1; touch make_spline
clean_spline:
	rm -f make_spline
	cd spline; make -j1 clean

make_kdtree:
	rm -f make_kdtree
	set -e; cd kdtree; make -j1; touch make_kdtree
clean_kdtree:
	rm -f make_kdtree
	cd kdtree; make -j1 clean

make_fson:
	rm -f make_fson
	set -e; cd fson; make -j1; touch make_fson
clean_fson:
	rm -f make_fson
	cd fson; make -j1 clean

make_slint:
	rm -f make_slint
	set -e; cd slint; make -j1; touch make_slint
clean_slint:
	rm -f make_slint
	cd slint; make -j1 clean

make_sondelib:
	rm -f make_sondelib
	set -e ; cd sondelib; make -j1; touch make_sondelib
clean_sondelib:
	rm -f make_sondelib
	cd sondelib ; make -j1 clean

make_slatec:
	rm -f make_slatec
	set -e; cd slatec; make -j1; touch make_slatec
clean_slatec:
	rm -f make_slatec
	cd slatec ; make -j1 clean

make_nr:
	rm -f make_nr
	set -e; cd nr; make -j1; touch make_nr
clean_nr:
	rm -f make_nr
	cd nr ; make -j1 clean

make_bacio:
	rm -f make_bacio
	set -e ; cd bacio ; make -j1 ; touch make_bacio
clean_bacio:
	rm -f make_bacio
	cd bacio ; make -j1 clean

make_g2c:
	rm -f make_g2c
	set -e ; cd g2c ; make -j1 ; touch make_g2c
clean_g2c:
	rm -f make_g2c
	cd g2c ; make -j1 clean

make_gctpc:
	rm -f make_gctpc
	set -e ; cd gctpc ; make -j1 ; touch make_gctpc
clean_gctpc:
	rm -f make_gctpc
	cd gctpc ; make -j1 clean

make_blas:
	rm -f make_blas
	set -e ; cd blas ; make -j1 ; touch make_blas
clean_blas:
	rm -f make_blas
	cd blas ; make -j1 clean

make_bufr: make_bacio make_sigio
	rm -f make_bufr
	set -e ; cd bufr ; make -j1 ; touch make_bufr
clean_bufr:
	rm -f make_bufr
	cd bufr ; make -j1 clean

make_sfcio:
	rm -f make_sfcio
	set -e ; cd sfcio ; make -j1 ; touch make_sfcio
clean_sfcio:
	rm -f make_sfcio
	cd sfcio ; make -j1 clean

make_sp: make_bacio
	rm -f make_sp
	set -e ; cd sp ; make -j1 ; touch make_sp
clean_sp:
	rm -f make_sp
	cd sp ; make -j1 clean

make_w3nco: make_w3emc
	rm -f make_w3nco
	set -e ; cd w3nco ; make -j1 ; touch make_w3nco
clean_w3nco:
	rm -f make_w3nco
	cd w3nco ; make -j1 clean

make_w3emc: make_sigio make_bacio
	rm -f make_w3emc
	set -e ; cd w3emc ; make -j1 ; touch make_w3emc
clean_w3emc:
	rm -f make_w3emc
	cd w3emc ; make -j1 clean

make_sigio:
	rm -f make_sigio
	set -e ; cd sigio ; make -j1 ; touch make_sigio
clean_sigio:
	rm -f make_sigio
	cd sigio ; make -j1 clean

make_hwrfutil:
	rm -f make_hwrfutil
	set -e ; cd hwrfutil ; make -j1 ; touch make_hwrfutil
clean_hwrfutil:
	rm -f make_hwrfutil
	cd hwrfutil ; make -j1 clean

make_g2:  make_bacio make_w3nco make_w3emc make_g2tmpl
	rm -f make_g2
	set -e ; cd g2 ; make -j1 ; touch make_g2
clean_g2:
	rm -f make_g2
	cd g2 ; make -j1 clean

make_g2tmpl:
	rm -f make_g2tmpl
	set -e ; cd g2tmpl ; make -j1 ; touch make_g2tmpl
clean_g2tmpl:
	rm -f make_g2tmpl
	cd g2tmpl ; make -j1 clean

make_ip: make_bacio
	rm -f make_ip
	set -e ; cd ip ; make -j1 ; touch make_ip
clean_ip:
	rm -f make_ip
	cd ip ; make -j1 clean

make_nemsio: 
	rm -f make_nemsio
	set -e ; cd nemsio ; make -j1 ; touch make_nemsio
clean_nemsio:
	rm -f make_nemsio
	cd nemsio ; make -j1 clean
