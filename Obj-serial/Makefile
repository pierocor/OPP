# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc 
CFLAGS=-Wall -pg -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include
#CFLAGS=-Wall -pg -g -no-pie -O3 -ffast-math -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/force.c ../src/input.c ../src/output.c ../src/utilities.c ../src/verlet1.c ../src/verlet2.c
VPATH=../src
OBJ=force.o input.o main.o output.o utilities.o verlet1.o verlet2.o

default: ../ljmd-serial.x

clean:
	rm -f *.mod *.o ../ljmd-serial.x

# linker rule
../ljmd-serial.x: $(OBJ) main.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<




main.o: ../main.c
	$(CC) -c $(CFLAGS) $<
