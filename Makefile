.PHONY: test clean

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

elphy: elphy.o driver.o io.o matrix.o supercell.o temperature.o sockets.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

%.dat: %.py
	python3 $<

test: elphy input.dat
	python3 test.py

clean:
	rm -f elphy *.o *.dat
