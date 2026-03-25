.PHONY: doc test ipi show clean clean_dat clean_ipi clean_all

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

model = graphene

elphy: elphy.o driver.o io.o matrix.o random.o sockets.o strain.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

doc: README.md
	pandoc -s -M pagetitle="elphy" -o index.html $<

test input.dat input.xyz: test.py elphy
	python3 $< $(model)

ipi ipi.pos_0.xyz: input.xml elphy input.dat input.xyz
	i-pi $< &
	sleep 3
	./elphy input.dat localhost:31415

show: ipi.pos_0.xyz
	python3 show.py

clean:
	rm -f elphy *.o

clean_dat:
	rm -f *.dat *.xyz

clean_ipi:
	rm -f RESTART \#* ipi.*

clean_all: clean clean_dat clean_ipi
