#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=18           # number of cores
#SBATCH --mem-per-cpu=5G             # shared memory, scaling with CPU request

module load miniconda
mamba activate snakemake
snakemake --cores $SLURM_CPUS_PER_TASK
