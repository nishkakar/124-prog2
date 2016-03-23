.PHONY: strassen
strassen:
	cc -std=c99 randmst.c -o randmst -lm
