include ../../configure.hwrf

OBJS	= hwrf_final_merge.o
EXEC	= ../../exec/hwrf_final_merge.exe

$(EXEC): $(OBJS)
	$(SFC) $(LDFLAGS) -o $@ $(OBJS)

.PHONY: clean

clean:
	$(RM) $(OBJS) *.exe

