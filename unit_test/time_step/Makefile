# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I $(HOME)/include
LDLIBS=-lm
HOME=../..

# list of source files
SRC=$(HOME)/src/*
VPATH=$(HOME)/src
OBJ=verlet1.o verlet2.o utilities.o single_time_step.o input.o

default: check

clean:
	rm -f *.mod *.o single_time_step.x single_time_step.dat

# linker rule
single_time_step.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

check: single_time_step.x
	./single_time_step.x
	cmp single_time_step.dat $(HOME)/reference/single_time_step.dat || exit 1
	rm single_time_step.dat single_time_step.x *.o
