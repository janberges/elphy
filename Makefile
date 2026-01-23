.PHONY: clean

elphy: elphy.c
	gcc -std=c89 -pedantic -Wall -o $@ $< -llapack -lblas -lm

clean:
	rm -f elphy
