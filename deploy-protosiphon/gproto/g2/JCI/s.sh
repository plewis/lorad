#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=jci
#SBATCH -o frt-jci.out
#SBATCH -e frt-jci.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

