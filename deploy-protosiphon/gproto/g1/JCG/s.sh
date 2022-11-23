#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=jcg
#SBATCH -o frt-jcg.out
#SBATCH -e frt-jcg.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

