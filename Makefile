.PHONY: test ipi show doc clean distclean

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall
LDLIBS = -llapack -lblas -lm

model = graphene

elphy: elphy.o driver.o io.o matrix.o random.o sockets.o strain.o supercell.o temperature.o
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}

%.o: %.c elphy.h
	$(CC) $(CFLAGS) -o $@ -c $<

test input.dat input.xyz: test.py elphy
	python3 $< $(model) input.dat input.xyz

ipi ipi.pos_0.xyz: input.xml elphy input.dat input.xyz
	i-pi $< &
	sleep 3
	./elphy input.dat localhost:31415

symmetric.xyz: elphy input.dat
	./$^ $@ 0

show: symmetric.xyz ipi.pos_0.xyz
	python3 show.py $^

doc:
	pandoc -s -M pagetitle="elphy" -o index.html README.md

clean:
	rm -f elphy *.o

distclean:
	rm -f $$(cat .gitignore)
