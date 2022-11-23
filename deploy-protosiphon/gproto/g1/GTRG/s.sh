#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=gtrg
#SBATCH -o frt-gtrg.out
#SBATCH -e frt-gtrg.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

