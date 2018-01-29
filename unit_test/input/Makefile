# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I $(HOME)/include
LDLIBS=-lm
HOME=../..

# list of source files
SRC=$(HOME)/src/input.c
VPATH=$(HOME)/src
OBJ=input.o input_param_test.o

default: check

clean:
	rm -f *.mod *.o input_param_test.x

# linker rule
input_param_test.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

check: input_param_test.x
	./input_param_test.x
	rm input_param_test.x *.o
