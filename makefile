CC = gcc
CFLAGS = -std=c99

all: hydroentropy hydroenthalpy

hydroentropy:
	$(CC) $(CFLAGS) src/hydroentropy/hydroentropy.c src/hydroentropy/hydrocluster.c src/hydroentropy/input.c src/hydroentropy/util.c src/hydroentropy/entropy.c src/hydroentropy/dbScan.c -lm -o bin/hydroentropy

hydroenthalpy:
	$(CC) $(CFLAGS) src/hydroenthalpy/hydroenthalpy.c src/hydroenthalpy/energy.c src/hydroenthalpy/energy_io.c -lm -o bin/hydroenthalpy 

clean:
	rm *.o 
