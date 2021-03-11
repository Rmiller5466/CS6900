#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj05S
#SBATCH --output=proj05S.txt

./proj05S.exe
