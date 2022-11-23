#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=gtr
#SBATCH -o frt-gtr.out
#SBATCH -e frt-gtr.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

