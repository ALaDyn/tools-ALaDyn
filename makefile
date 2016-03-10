
TOOLS := analizza_binario \
leggi_bin


all:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i ;\
	done

clean:
	@for i in $(TOOLS) ; do \
	$(MAKE) -C $$i clean ;\
	done

