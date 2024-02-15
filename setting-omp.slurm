#!/bin/bash
#SBATCH --job-name=floyd-warshall
#SBATCH --output=out/test-omp-32.out
#SBATCH --error=err/test-omp.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000
#SBATCH --cpu-freq=performance
#SBATCH --partition=skyvolta

module load cuda/10.2

make type=omp

ulimit -s $((10*1024*1024))

export OMP_CANCELLATION=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun ./heat2D-omp
