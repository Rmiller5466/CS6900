#!/usr/bin/bash

# compile the source
/act/openmpi-2.0/gcc/bin/mpicc -o $1.omp $1.c -lm -fopenmp

# make it executable
chmod 700 $1.omp

