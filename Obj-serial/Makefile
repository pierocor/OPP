# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/force.c ../src/input.c ../src/output.c ../src/utilities.c ../src/verlet1.c ../src/verlet2.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x

clean:
	rm -f *.mod *.o ../ljmd-serial.x

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ main.o $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

main.o: $(OBJ)

../ljmd-serial.x: main.o

main.o: ../main.c
	$(CC) -c $(CFLAGS) $<
