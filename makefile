
TOOLS := analizza_binario \
combina_vtk \
exponential_fit \
fix_nptot \
InGenUO \
interpolate_scan_results \
leggi_bin \
leggi_diagspec \
leggi_partdist \
lightweight_coredump_analyzer \
logaritmic_fit


all:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i ;\
	done

clean:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i clean ;\
	done

cleanall: clean
	