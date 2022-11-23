#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3jcig
#SBATCH -o frt-3jcig.out
#SBATCH -e frt-3jcig.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

