#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=gtri
#SBATCH -o frt-gtri.out
#SBATCH -e frt-gtri.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
$HOME/yubo/loradmth

