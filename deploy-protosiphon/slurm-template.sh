#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=__PREFIX__
#SBATCH -o frt-__PREFIX__.out
#SBATCH -e frt-__PREFIX__.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

