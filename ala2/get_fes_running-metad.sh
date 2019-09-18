#! /bin/bash

kbt=2.494339
bck=''
#bck='bck.0.'
stride=20

tot=`ls Hills.*.data |wc -l`
for i in  `seq 0 $((tot-1))`
do
  echo -en "  summing... \r"
  plumed sum_hills --mintozero --hills ${bck}Hills.$i.data --min -pi,-pi --max pi,pi --bin 100,100 --stride $stride > log.sum_hills
  bck.meup.sh fes_running.$i fes_running.${i}.data
  mkdir fes_running.$i
  mv fes_*.dat fes_running.$i
  mv `ls -v fes_running.$i/* |tail -1` fes_running.${i}.data

  bck.meup.sh fes_deltaF.${i}.data
  echo "#time deltaF" > fes_deltaF.${i}.data
  t=0
  for f in `ls -v fes_running.${i}/*`
  do
    t=$((t+stride))
    echo -en "  working... t=$t \r"
    awk -v kbt=$kbt -v time=$t '{if($1!="#!" && NF>1) {if($1>0 && $1<2.3) basinB+=exp(-$3/kbt); else basinA+=exp(-$3/kbt);} }END{print time,kbt*log(basinA/basinB)}' $f >> fes_deltaF.${i}.data
    mv $f fes_running.${i}/fes_running.t-${t}.data
  done
done

outfile=stats-fes_deltaF.data
bck.meup.sh -i $outfile

paste ${bck}fes_deltaF.*.data |\
  awk 'BEGIN{print "#average std_dev"}
       NR>1{
         av=0; av2=0; 
         for (i=2; i<=NF; i+=2) {av+=$i; av2+=$i^2;};
         av*=2./NF; av2*=2./NF; 
         print $1,av,sqrt(av2-av^2)
       }' > $outfile
