#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj05
#SBATCH --output=proj05.txt

# cs6900
# Project 05
# Ryan Miller
# w051rem
# Due 19 March 2021
# System = bender
# Job Control File = proj05.sbatch
# Results file     = proj05.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi-2.0/gcc

#mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 2 ./proj05.mpi -r 1600 -R 10
#mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 3 ./proj05.mpi -f testFile.txt -F testref.txt
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj05.mpi -r 500 -R 10
echo 
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj05.mpi -r 1000 -R 10
echo 
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj05.mpi -r 2500 -R 10
echo 
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 16 ./proj05.mpi -r 5000 -R 10

