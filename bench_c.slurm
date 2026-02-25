#!/bin/sh

#SBATCH --time=3:00:0
#SBATCH --ntasks=1
#SBATCH --partition=serial
#SBATCH --qos=serial

module load PrgEnv-gnu

# Switch to current working directory
cd C

# Run the parallel program
./MD



