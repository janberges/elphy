.PHONY: test ipi show clean clean_dat clean_ipi clean_all

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

elphy: elphy.o driver.o io.o matrix.o supercell.o temperature.o sockets.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

input.dat: data.py
	python3 $<

test: elphy input.dat
	python3 test.py

ipi ipi.pos_0.xyz: elphy input.dat
	./elphy input.dat input.xyz
	i-pi input.xml &
	sleep 2
	./elphy input.dat

show: ipi.pos_0.xyz
	python3 show.py

clean:
	rm -f elphy *.o

clean_dat:
	rm -f *.dat *.xyz

clean_ipi:
	rm -f RESTART \#* ipi.*

clean_all: clean clean_dat clean_ipi
