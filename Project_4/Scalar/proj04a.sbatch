#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj04a
#SBATCH --output=proj04a.txt

# cs6900
# Project 04a
# Ryan Miller
# w051rem
# Due 1 March 2021
# System = bender
# Job Control File = proj04a.sbatch
# Results file     = proj04a.txt

./proj04a.exe 20 20 1 -1 135 135

