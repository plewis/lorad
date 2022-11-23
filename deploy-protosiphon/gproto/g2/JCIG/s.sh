#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=jcig
#SBATCH -o frt-jcig.out
#SBATCH -e frt-jcig.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

