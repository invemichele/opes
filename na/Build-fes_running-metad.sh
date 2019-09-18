#!/bin/bash

n_walkers=`ls Colvar.*data |wc -l`
stride_time=2000
echo " n_walkers = $n_walkers"

plumed sum_hills --mintozero --hills Hills.data --min 1.5 --max 2.7 --bin 99 --stride $((stride_time*n_walkers)) > log.sum_hills
bck.meup.sh fes_running
mkdir fes_running
for f in `ls -v fes_*.dat`
do
  n=${f#fes_};n=${n%.dat};n=$((n+1))
  echo -en " $f\r"
  mv $f fes_running/fes_running.t-$((n*stride_time)).data
done
mv fes_running/fes_running.t-$((n*stride_time)).data fes_running.data
