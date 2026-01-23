.PHONY: test clean

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

elphy: elphy.o io.o matrix.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: elphy
	python3 test.py

clean:
	rm -f elphy *.o
