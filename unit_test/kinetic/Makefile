# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I $(HOME)/include
LDLIBS=-lm
HOME=../..

# list of source files
SRC=$(HOME)/src/utilities.c
VPATH=$(HOME)/src
OBJ=kinetic_test.o utilities.o

default: check

clean:
	rm -f *.mod *.o kinetic_test.x kinetic_test.dat

# linker rule
kinetic_test.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

check: kinetic_test.x
	./kinetic_test.x
	cmp kinetic_test.dat $(HOME)/reference/kinetic_test.dat || exit 1
	rm kinetic_test.dat kinetic_test.x *.o
