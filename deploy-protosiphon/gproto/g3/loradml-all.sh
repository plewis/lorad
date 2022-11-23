#!/bin/bash

for MODEL in JC JCI JCG JCIG GTR GTRI GTRG GTRIG 3JC 3JCI 3JCG 3JCIG 3GTR 3GTRI 3GTRG 3GTRIG
do
    cd $MODEL
    LCMODEL="$(tr [A-Z] [a-z] <<< "$MODEL")"
    python3 ../../filter.py $LCMODEL-logtransformed-params.txt $LCMODEL-logtransformed-params.txt 
    $HOME/yubo/loradML > ../loradML-$MODEL-output.txt
    cd ..
done
wait
