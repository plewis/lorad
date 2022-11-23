#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=3gtrig
#SBATCH -o frt-3gtrig.out
#SBATCH -e frt-3gtrig.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

