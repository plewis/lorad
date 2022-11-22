#!/bin/bash

mkdir g

INDEXES=({1..20})
for i in ${INDEXES[@]} ; do
    python3 deploy.py $i $i`pickseed 3`
    mv g$i g
done

cp go-example.sh g/go.sh
cp coverage-series-example.sh g/coverage-series.sh
