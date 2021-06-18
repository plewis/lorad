#!/bin/bash

#SBATCH --job-name=hpdml
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=paul.lewis@uconn.edu
#SBATCH -o marglike%j.out
#SBATCH -e marglike%j.err

LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
./hpdml
