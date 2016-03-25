.PHONY: strassen
strassen:
	gcc -std=c99 strassen3.c -o strassen -lm
