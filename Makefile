.PHONY: test clean

elphy: elphy.c
	gcc -std=c89 -pedantic -Wall -o $@ $< -llapack -lblas -lm

test: elphy
	python3 test.py

clean:
	rm -f elphy
