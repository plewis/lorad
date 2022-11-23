#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3jcg
#SBATCH -o frt-3jcg.out
#SBATCH -e frt-3jcg.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

