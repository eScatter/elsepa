FC = f77
FCFLAGS = -O3

PROGRAMS = elscata elscatm

all: $(PROGRAMS)

elscata.o: src/getpath.f

elscatm.o: src/getpath.f

elscata: elscata.o elsepa.o

elscatm: elscatm.o elsepa.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: src/%.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o

#PREFIX is environment variable, but if it is not set, then set default value
PREFIX ?= /usr/local

.PHONY: install
install: all
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $(PROGRAMS) $(DESTDIR)$(PREFIX)/bin
	mkdir -p $(DESTDIR)$(PREFIX)/share/elsepa
	cp -r data/ $(DESTDIR)$(PREFIX)/share/elsepa
	cp -r examples/ $(DESTDIR)$(PREFIX)/share/elsepa

	@echo ""
	@echo "*******************************************************************************"
	@echo "  BEFORE RUNNING ELSEPA, PLEASE SET THE ELSEPA_DATA ENVIRONMENT VARIABLE"
	@echo "  TO "$(DESTDIR)$(PREFIX)/share/elsepa/data"."
	@echo ""
	@echo "  For example, put the following line in your .bashrc:"
	@echo "  export ELSEPA_DATA=$(DESTDIR)$(PREFIX)/share/elsepa/data"
	@echo "*******************************************************************************"

.PHONY: uninstall
uninstall:
	for p in $(PROGRAMS) ; do rm -f $(DESTDIR)$(PREFIX)/bin/$$p ; done
	rm -rf $(DESTDIR)$(PREFIX)/share/elsepa
