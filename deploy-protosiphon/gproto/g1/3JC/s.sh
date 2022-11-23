#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3jc
#SBATCH -o frt-3jc.out
#SBATCH -e frt-3jc.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

