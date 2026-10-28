.PHONY: test ipi ipi_unix ipi_shm md metad show_ipi show_md show_metad clean distclean lint

CC := gcc
CFLAGS := -ansi -Wpedantic -Wall -Wextra
LDFLAGS := $(CFLAGS)
LDLIBS := -llapack -lblas -lm

model := graphene
units := Ha

elphy: elphy.o driver.o io.o matrix.o random.o strain.o supercell.o temperature.o sockets.o plumed.o
	$(strip $(CC) $(LDFLAGS)) -o $@ $^ $(strip $(LDLIBS) $(EXTRA_LDLIBS))
	test -f $@.exe -a ! -f $@ && printf '#!/bin/sh\nwine ./$@.exe $$@' > $@ && chmod +x $@

%.o: %.c elphy.h
	$(strip $(CC) $(CFLAGS)) -o $@ -c $<

test input.dat mdargs.dat plumed.dat: test.py elphy
	python3 $< $(model) $(units) input.dat mdargs.dat plumed.dat

input.xyz: elphy input.dat
	./$^ -1 0.1 > $@

symmetric.xyz: elphy input.dat
	./$^ -1 0 > $@

ipi ipi.pos_0.xyz: input.xml elphy input.dat input.xyz
	i-pi $< & ./elphy input.dat localhost:31415

ipi_unix: input.xml elphy input.dat input.xyz
	sed s/inet/unix/ $< > input_unix.xml
	i-pi input_unix.xml & ./elphy input.dat localhost

ipi_shm: input.xml elphy input.dat input.xyz
	sed s/inet/shm/ $< > input_shm.xml
	i-pi input_shm.xml & ./elphy input.dat localhost/shm

md md.xyz nvt.xyz nve.xyz dmp.xyz: elphy input.dat
	./$^ 1000:10 $$(cat mdargs.dat) < /dev/null > md.xyz 2> nvt.xyz
	./$^ 1000:10 $$(awk '$$2=0;1' mdargs.dat) < nvt.xyz >> md.xyz 2> nve.xyz
	./$^ 1000:10 $$(awk '$$3=0;1' mdargs.dat) < nve.xyz >> md.xyz 2> dmp.xyz

metad metad.xyz plumed.log: elphy input.dat plumed.dat
	./elphy input.dat 5000:10 $$(cat mdargs.dat) plumed.dat:plumed.log < /dev/null > metad.xyz

show_ipi: show.py symmetric.xyz ipi.pos_0.xyz
	python3 $^

show_md: show.py symmetric.xyz md.xyz
	python3 $^

show_metad: show.py symmetric.xyz metad.xyz
	python3 $^

index.html: README.md
	pandoc -s -M pagetitle="elphy" -o $@ $<

%.png: %.svg
	inkscape -w 600 -o $@ $<
	python -c "import storylines as sl; sl.save('$@', sl.load('$@'))"

clean:
	rm -f elphy elphy.exe *.o

distclean:
	rm -f $$(cat .gitignore)

lint:
	clang-tidy *.[ch] -- $(CFLAGS)
