.PHONY: strassen
strassen:
	cc -std=c99 strassen.c -o strassen -lm
