#!/bin/bash

INDEXES=({1..3})
for i in ${INDEXES[@]} ; do
  cd g3
  . submit-all.sh
  cd ..
done
