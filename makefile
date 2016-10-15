
TOOLS := binary_analyzer \
binary_decoder \
diagspec_reader \
exponential_fit \
extract_nptot \
fix_nptot \
InGenUO \
interpolate_scan_results \
lightweight_coredump_analyzer \
logaritmic_fit \
merge_vtk \
partdist_reader


all:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i ;\
	done

clean:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i clean ;\
	done

cleanall: clean
	
