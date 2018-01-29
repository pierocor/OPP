# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fPIC -fomit-frame-pointer -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/
VPATH=../src
OBJ=force.o input.o output.o utilities.o verlet1.o verlet2.o

default: shared

clean:
	rm -f *.so *.mod *.o
shared: py_shared.so

py_shared.so: $(OBJ)
	$(CC) -shared $(OBJ) -o py_shared.so

%.o: %.c
	$(CC) -c $(CFLAGS) $<
