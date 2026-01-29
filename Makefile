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

sockets.c:
	wget https://raw.githubusercontent.com/i-pi/i-pi/refs/heads/main/drivers/f90/$@

sockets.o: sockets.c
	$(CC) -o $@ -c $<

test: elphy input.dat
	python3 test.py

clean:
	rm -f elphy *.o *.dat sockets.c
