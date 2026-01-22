.PHONY: clean

elphy: elphy.c
	gcc -std=c89 -pedantic -Wall -o $@ $< -llapack -lblas

clean:
	rm -f elphy
