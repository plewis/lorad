#!/bin/bash

#SBATCH --partition=SkylakePriority
#SBATCH --job-name=covbycodon
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --output=bycodon_%A-%a.out  # Standard output and error log
#SBATCH --array=1-11                # Array range

# Copy this script and change "bycodon" everywhere to "unpart", "bygene", 
# or "byboth" as desired. Also ensure loradML executable is in $PATH or 
# change last line to provide fully-qualified path to the loradML executable.

COVERAGES=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99)

cd g1/bycodon/lorad
cp loradml-untransformed.conf loradml.conf
phi=${COVERAGES[($SLURM_ARRAY_TASK_ID - 1)]}
echo "task $SLURM_ARRAY_TASK_ID: $i (coverage = $phi)"
export TIMEFORMAT="user-seconds %3U"
time loradML --coverage=$phi --mcse --tbratio=10
