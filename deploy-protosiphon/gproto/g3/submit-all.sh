#!/bin/bash

for MODEL in JC JCI JCG JCIG GTR GTRI GTRG GTRIG 3JC 3JCI 3JCG 3JCIG 3GTR 3GTRI 3GTRG 3GTRIG
do
    cd $MODEL
    sbatch s.sh
    cd ..
done
wait
