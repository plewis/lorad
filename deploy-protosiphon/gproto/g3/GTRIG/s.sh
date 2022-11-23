#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=gtrig
#SBATCH -o frt-gtrig.out
#SBATCH -e frt-gtrig.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

