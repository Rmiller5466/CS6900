#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj04a
#SBATCH --output=proj04a.txt
echo "1 1 10 80"
./proj04a.exe 100 100 1 1 0 45
echo
#echo "2 -2 10 80"
#./proj04a.exe 10 10 2 -2 10 80
#echo 
#echo "2 -2 45 45"
#./proj04a.exe 10 10 2 -2 45 45
#echo 
#echo "2 2 45 45"
#./proj04a.exe 10 10 2 2 45 45
#echo 
