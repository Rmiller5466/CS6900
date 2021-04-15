#!/usr/bin/bash


if [ $# \> 1 ];
then
    # compile the source
    /act/openmpi-2.0/gcc/bin/mpicc -o $1.mpi $1.c $2.c -lm -fopenmp
else
    /act/openmpi-2.0/gcc/bin/mpicc -o $1.mpi $1.c -lm -fopenmp
fi
# make it executable
chmod 700 $1.mpi

