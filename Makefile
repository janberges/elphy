.PHONY: test ipi show doc clean distclean

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

symmetric.xyz: elphy input.dat
	./$^ -1 0 > $@

show: symmetric.xyz ipi.pos_0.xyz
	python3 show.py $^

doc:
	pandoc -s -M pagetitle="elphy" -o index.html README.md

%.png: %.svg
	inkscape -w 600 -o $@ $<
	python -c "import storylines as sl; sl.save('$@', sl.load('$@'))"

clean:
	rm -f elphy *.o

distclean:
	rm -f $$(cat .gitignore)
