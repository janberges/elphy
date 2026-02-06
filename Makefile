.PHONY: test ipi ipi_elphmod show clean clean_dat clean_ipi clean_all

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

elphy: elphy.o driver.o io.o matrix.o random.o sockets.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

input.dat driver.pickle: data.py
	python3 $<

test: elphy input.dat driver.pickle
	python3 test.py

input.xyz: elphy input.dat
	./$^ $@ 0.1

ipi ipi.pos_0.xyz: input.xml elphy input.dat input.xyz
	i-pi $< &
	sleep 3
	./elphy input.dat localhost:31415

ipi_elphmod: input.xml input.xyz driver.pickle
	i-pi $< &
	sleep 3
	i-pi-driver-py -m elphmod -o driver=driver.pickle -a localhost -p 31415

show: ipi.pos_0.xyz driver.pickle
	python3 show.py

clean:
	rm -f elphy *.o

clean_dat:
	rm -f *.dat *.xyz

clean_ipi:
	rm -f RESTART \#* ipi.*

clean_all: clean clean_dat clean_ipi
