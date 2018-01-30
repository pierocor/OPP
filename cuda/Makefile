# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include
#CFLAGS=-Wall -g -no-pie -O3 -ffast-math -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/force.cu ../src/input.c ../src/output.c ../src/utilities.c ../src/verlet1.c ../src/verlet2.c
VPATH=../src
OBJ=forces.o input.o main.o output.o utilities.o verlet1.o verlet2.o

default: ./cuda.x

clean:
	rm -f *.mod *.o

# linker rule
cuda.x: $(OBJ) main.o
	nvcc -c ../src/force.cu forces.o -I ../include/

# compilation pattern rule for objects
%.o: %.c
	nvcc -c $(CFLAGS) $<

forces.o: ../src/force.cu
	nvcc -c ../src/force.cu forces.o -I ../include/

main.o: ./cuda_main.cu
	nvcc -c ../src/force.cu forces.o -I ../include/
