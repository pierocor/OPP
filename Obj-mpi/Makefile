# -*- Makefile -*-
SHELL=/bin/sh
CC=mpicc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -DMPI -I ../include
LDLIBS=-lm

# list of source files
SRC=force.c input.c output.c utilities.c verlet1.c verlet2.c
VPATH=../src
OBJ=force.o input.o main.o output.o utilities.o verlet1.o verlet2.o

default: ../ljmd-mpi.x

clean:
	rm -f *.mod *.o ../ljmd-mpi.x

# linker rule
../ljmd-mpi.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

main.o: ../main.c
	$(CC) -c $(CFLAGS) $< -o $@
