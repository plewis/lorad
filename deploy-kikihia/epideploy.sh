#!/bin/bash

mkdir g

INDEXES=({1..100})
for i in ${INDEXES[@]} ; do
    python3 deploy.py $i $i`pickseed 3`
    mv g$i g
done

cp go-template.sh g/go.sh
