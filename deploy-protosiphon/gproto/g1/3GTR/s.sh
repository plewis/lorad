#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3gtr
#SBATCH -o frt-3gtr.out
#SBATCH -e frt-3gtr.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

