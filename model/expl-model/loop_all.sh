#!/bin/bash

for i in `seq 0 9`
do
  dir=rep_$i
  mkdir $dir
  cd $dir
  cp ../../inputs/md_input.${i}.dat md_input.dat
  cp ../../inputs/md_potential.dat ../plumed.dat .
  ../queue_md.sh
  cd ..
done
