#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj02
#SBATCH --output=array%A_%a.txt
#SBATCH --array=1-16

for i in {1..512}; do ./decrypt ${SLURM_ARRAY_TASK_ID} $i w051rem; done
