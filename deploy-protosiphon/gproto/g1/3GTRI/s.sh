#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3gtri
#SBATCH -o frt-3gtri.out
#SBATCH -e frt-3gtri.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

