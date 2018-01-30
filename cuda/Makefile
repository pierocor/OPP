# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -I ../include
#CFLAGS=-Wall -g -no-pie -O3 -ffast-math -std=c11 -I ../include
LDLIBS=-lm

# list of source files
SRC=../src/input.c ../src/output.c ../src/utilities.c ../src/verlet1.c ../src/verlet2.c
VPATH=../src
OBJ=cuda.o input.o cuda_main.o output.o utilities.o verlet1.o verlet2.o

default: ./cuda.x

clean:
	rm -f *.mod *.o

# linker rule
cuda.x: $(OBJ) cuda_main.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

cuda.o: ./cuda.cu
	nvcc -c ./cuda.cu -I ../include/

cuda_main.o: ./cuda_main.cu
	nvcc -c cuda_main.cu -I ../include/
