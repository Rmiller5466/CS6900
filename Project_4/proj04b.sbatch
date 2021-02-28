#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj4b
#SBATCH --output=proj4b.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi-2.0/gcc

echo
echo "Arguments: 10^9 0 1 -1 45 45" 
echo
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj04b.mpi 1000000000 0 1 -1 45 45
