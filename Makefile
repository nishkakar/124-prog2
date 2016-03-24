.PHONY: strassen
strassen:
	cc -std=c99 strassen2.c -o strassen -lm
