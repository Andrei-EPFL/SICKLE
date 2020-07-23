#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J G2048CIC 	# Name of the job
#SBATCH -p p4		# Partition
#SBATCH -N 1            # number of nodes
#SBATCH -c 16           # number of cpus per tasks
#SBATCH -o ./output/out_linhaloCIC.out
#SBATCH -e ./output/err_linhaloCIC.err

export OMP_NUM_THREADS=16

srun ./linhalo -c ./linhalo.conf
