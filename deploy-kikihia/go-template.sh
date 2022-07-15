#!/bin/bash

INDEXES=(64 65 66 67 68 69 70 71 72 73 74 75)
for i in ${INDEXES[@]} ; do
    cd g$i
    . submit-ghme.sh
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