#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3gtrg
#SBATCH -o frt-3gtrg.out
#SBATCH -e frt-3gtrg.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

