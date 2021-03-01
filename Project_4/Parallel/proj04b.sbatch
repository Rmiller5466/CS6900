#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj04b
#SBATCH --output=proj04b.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi-2.0/gcc

mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj04b.mpi
