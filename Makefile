.PHONY: strassen
strassen:
	cc -std=c99 strassen3.c -o strassen -lm
