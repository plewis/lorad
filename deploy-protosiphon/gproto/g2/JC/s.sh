#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=jc
#SBATCH -o frt-jc.out
#SBATCH -e frt-jc.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

