# -*- Makefile -*-
SHELL=/bin/sh
CC= gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=force.c input.c output.c utilities.c verlet1.c verlet2.c
VPATH=../src
OBJ=force.o input.o main.o output.o utilities.o verlet1.o verlet2.o

default: mpi

mpi: CFLAGS = -Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include -DMPI
mpi: CC = mpicc
mpi: ../ljmd-mpi.x

omp: CFLAGS = -Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include -DOPENMP -fopenmp
omp: ../ljmd-omp.x

combo: CFLAGS = -Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include -DMPI -DOPENMP -fopenmp
combo: CC = mpicc
combo: ../ljmd-combo.x

clean:
	rm -f *.mod *.o ../ljmd-mpi.x ../ljmd-omp.x ../ljmd-combo.x

# linker rule

../ljmd-mpi.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)


../ljmd-omp.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

../ljmd-combo.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)
# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

main.o: ../main.c
	$(CC) -c $(CFLAGS) $< -o $@
