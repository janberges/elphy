.PHONY: test ipi ipi_unix ipi_shm show_ipi md show_md clean distclean

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall -Wno-parentheses
LDLIBS = -llapack -lblas -lm

model = graphene
units = Ha

elphy: elphy.o driver.o io.o matrix.o random.o sockets.o strain.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

test input.dat: test.py elphy
	python3 $< input.dat $(model) $(units)

input.xyz: elphy input.dat
	./$^ -1 0.1 > $@

ipi ipi.pos_0.xyz: input.xml elphy input.dat input.xyz
	i-pi $< &
	sleep 3
	./elphy input.dat localhost:31415

ipi_unix: input.xml elphy input.dat input.xyz
	sed s/inet/unix/ $< > input_unix.xml
	i-pi input_unix.xml &
	sleep 3
	./elphy input.dat localhost

ipi_shm: input.xml elphy input.dat input.xyz
	sed s/inet/shm/ $< > input_shm.xml
	i-pi input_shm.xml &
	sleep 3
	./elphy input.dat localhost/shm

symmetric.xyz: elphy input.dat
	./$^ -1 0 > $@

show_ipi: symmetric.xyz ipi.pos_0.xyz
	python3 show.py $^

nvt = $$(tail -n 3 input.dat | head -n 1) < /dev/null > md.xyz 2> nvt.xyz
nve = $$(tail -n 2 input.dat | head -n 1) < nvt.xyz >> md.xyz 2> nve.xyz
dmp = $$(tail -n 1 input.dat) < nve.xyz >> md.xyz 2> dmp.xyz

md md.xyz nvt.xyz nve.xyz dmp.xyz: elphy input.dat
	@echo "$(nvt)"; $(nvt)
	@echo "$(nve)"; $(nve)
	@echo "$(dmp)"; $(dmp)

show_md: symmetric.xyz md.xyz
	python3 show.py $^

index.html: README.md
	pandoc -s -M pagetitle="elphy" -o $@ $<

%.png: %.svg
	inkscape -w 600 -o $@ $<
	python -c "import storylines as sl; sl.save('$@', sl.load('$@'))"

clean:
	rm -f elphy *.o

distclean:
	rm -f $$(cat .gitignore)
