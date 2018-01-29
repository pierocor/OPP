# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I $(HOME)/include
LDLIBS=-lm
HOME=../..

# list of source files
SRC=$(HOME)/src/force.c $(HOME)/src/utilities.c
VPATH=$(HOME)/src
OBJ=force.o force_test.o utilities.o

default: check

clean:
	rm -f *.mod *.o force_test.x force_test.dat

# linker rule
force_test.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

check: force_test.x
	./force_test.x
	cmp force_test.dat $(HOME)/reference/force_test.dat || exit 1
	rm force_test.dat force_test.x *.o
