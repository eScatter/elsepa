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

