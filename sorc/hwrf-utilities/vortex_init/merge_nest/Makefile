include ../../configure.hwrf
include ../../pure-openmp.inc

OBJS	= merge_nest_4x_step12_3n.o merge_nest_4x_step12_enkf.o interp_coef.o

ALL: hwrf_merge_nest_4x_step12_3n.exe hwrf_merge_nest_4x_step12_enkf.exe

.PHONY: ALL

hwrf_merge_nest_4x_step12_3n.exe: merge_nest_4x_step12_3n.o interp_coef.o
	$(SFC) $(LDFLAGS) -o $@ merge_nest_4x_step12_3n.o interp_coef.o
	$(CP) *.exe ../../exec
	$(RM) *.exe

hwrf_merge_nest_4x_step12_enkf.exe: merge_nest_4x_step12_enkf.o 
	$(SFC) $(LDFLAGS) -o $@ merge_nest_4x_step12_enkf.o
	$(CP) *.exe ../../exec
	$(RM) *.exe

.PHONY: clean

clean:
	$(RM) $(OBJS) *.exe

