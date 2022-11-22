#!/bin/bash

# To submit 20 runs stored in directories named g$i
# Uncomment appropriate submit-xxx.sh line
INDEXES=({1..20})
for i in ${INDEXES[@]} ; do
    cd g$i
    . submit-lorad.sh
    #. submit-gss.sh
    cd ..
done

# To run loradML on parameter files in 20 directories named g$i/unpart/lorad
# Change unpart to bygene, bycodon, or byboth as needed
#INDEXES=({1..20})
#for i in ${INDEXES[@]} ; do
#    cd g$i/unpart/lorad
#    cp loradml-untransformed.conf loradml.conf
#    /home/pol02003/yubo/loradML --ghm --quiet --coverage=0.5
#    cd ../../..
#done

# To harvest results from 20 generalized steppingstone runs
# Change unpart to bygene, bycodon, or byboth as needed
#INDEXES=({1..20})
#for i in ${INDEXES[@]} ; do
#    cd g$i/unpart/gss
#    egrep "log\(marginal likelihood\) =" ./* --include="*.out"
#    cd ..
#done
#