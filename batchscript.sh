#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J G2048CIC 	# Name of the job
#SBATCH -p p5		# Partition
#SBATCH -N 1            # number of nodes
#SBATCH -c 32           # number of cpus per tasks
#SBATCH -o ./output/out.linhaloCIC.out
#SBATCH -e ./output/err.linhaloCIC.err

export OMP_NUM_THREADS=32

srun ./linhalo -c ./linhalo.conf
