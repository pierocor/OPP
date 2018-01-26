# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-shared -fPIC -Wall -O -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/
VPATH=../src
OBJ=force.so input.so output.so utilities.so verlet1.so verlet2.so

clean:
	rm -f *.so *.mod *.o
shared: $(OBJ)

%.so: ../src/%.c
	$(CC) -o $@ $(CFLAGS) $<

#gcc -o hello.so -shared -fPIC hello.c
