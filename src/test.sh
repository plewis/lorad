#!/bin/bash

module load gcc/10.2.0
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
./hpdml --help
