#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj04a
#SBATCH --output=proj04a.txt

./proj04a.exe 10 10 2 3 20 80
echo
#./proj04a.exe 6 6 2 -1 60 45
#echo 
#./proj04a.exe 6 9 1 -1 45 45
#echo 
