#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=4
#SBATCH --job-name=proj03
#SBATCH --output=proj03.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi-2.0/gcc

# Run on a single Node
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 4 ./proj03.mpi 1000
