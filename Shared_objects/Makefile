# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fPIC -fomit-frame-pointer -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/
VPATH=../src
OBJ=verlet1.o force.o input.o output.o utilities.o  verlet2.o

default: shared

clean:
	rm -f *.so *.mod *.o *.dat *.xyz
shared: py_shared.so

py_shared.so: $(OBJ)
	$(CC) -shared $(OBJ) -o py_shared.so

%.o: %.c
	$(CC) -c $(CFLAGS) $<
