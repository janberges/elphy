.PHONY: test clean

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

elphy: elphy.o io.o matrix.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

el.dat ph.dat elph.dat: data.py
	python3 $<

test: elphy el.dat ph.dat elph.dat
	python3 test.py

clean:
	rm -f elphy *.o *.dat
