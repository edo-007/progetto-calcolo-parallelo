#!/bin/bash
#SBATCH --job-name=floyd-warshall
#SBATCH --output=out/test-serial.out
#SBATCH --error=err/test-serial.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --cpu-freq=performance
#SBATCH --partition=skyvolta

module load devtoolset-9

make type=serial

ulimit -s $((10*1024*1024))

srun ./heat2D-serial
