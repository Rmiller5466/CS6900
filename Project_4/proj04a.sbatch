#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj04a
#SBATCH --output=proj04a.txt

./proj04a.exe 6 6 3 -1 45 45
