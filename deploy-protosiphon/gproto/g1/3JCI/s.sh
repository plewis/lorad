#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3jci
#SBATCH -o frt-3jci.out
#SBATCH -e frt-3jci.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

