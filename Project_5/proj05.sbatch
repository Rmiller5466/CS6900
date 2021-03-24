#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=5
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

mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 3 ./proj05.mpi -d 5 -F testref.txt
mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 3 ./proj05.mpi -f testFile.txt -F testref.txt
#mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 3 ./proj05.mpi -d 10 -R 3
#echo 
#mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 5 ./proj05.mpi -r 15 -R 3
#mpirun --mca btl_tcp_if_include en0 --mca btl tcp,self -np 5 ./proj05.mpi -r 7 -R 3

