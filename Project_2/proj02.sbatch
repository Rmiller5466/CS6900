#!/bin/bash
#
# cs6900-01
# Project 02
# Ryan Miller
# w051rem
# Due 03 Feb 2021
# System = bender
# Job Control File = proj02.sbatch
# Results File = proj02.txt
#
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj02
#SBATCH --output=array%A_%a.txt
#SBATCH --array=1-16

# Simple for loop that will test all 512 combinations with the given ID(1-16)
for i in {1..512}; do ./decrypt ${SLURM_ARRAY_TASK_ID} $i w051rem; done