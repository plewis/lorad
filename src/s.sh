#!/bin/bash

#SBATCH --job-name=hpdml
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=analisa.milkey@uconn.edu
#SBATCH -o marglike%j.out
#SBATCH -e marglike%j.err

module load gcc/10.2.0
LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
./hpdml
