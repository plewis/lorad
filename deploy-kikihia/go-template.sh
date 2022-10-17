#!/bin/bash

INDEXES=({1..20})
for i in ${INDEXES[@]} ; do
    cd g$i
    . submit-lorad.sh
    cd ..
done

#mkdir s
#for i in ${INDEXES[@]} ; do
#    cd g$i
#    python summary.py
#    cp output-summary.txt ../s/g$i.txt
#    cd ..
#done
#