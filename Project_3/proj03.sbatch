#!/bin/bash

#SBATCH -N 4
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj3TestRun
#SBATCH --output=proj3Test.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi-2.0/gcc

# Run on a single Node
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 4 ./proj03.mpi 100000 >> results.txt
