#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J pkvoid 	# Name of the job
#SBATCH -p p4		# Partition
#SBATCH -N 1            # number of nodes
#SBATCH -c 16           # number of cpus per tasks
#SBATCH -o ./output/out.1dive.out
#SBATCH -e ./output/err.1dive.err

export OMP_NUM_THREADS=16

srun temp_DIVE.txt
