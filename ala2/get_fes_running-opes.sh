#! /bin/bash

#kbt=2.494339
bck=''
#bck='bck.0.'
stride=20

tot=`ls Kernels.*.data |wc -l`
for i in  `seq 0 $((tot-1))`
do
  echo " getting prob ${i}..."
  bck.meup.sh Prob_running.${i}
  mkdir Prob_running.${i}
  awk -v i=$i '{if ($2=="FIELDS") outname=sprintf("Prob_running.%d/prob_running-%d.data",i,++n); print $0 > outname }' ${bck}Prob.${i}.data

  bck.meup.sh fes_running.${i} fes_running.${i}.data fes_deltaF.${i}.data
  echo "#time deltaF" > fes_deltaF.${i}.data
  mkdir fes_running.${i}
  t=0
  for f in `ls -v Prob_running.${i}/*`
  do
    t=$((t+stride))
    outfile="fes_running.${i}/fes_running.t-${t}.data"
    echo -en " $outfile \r"
    ./Build-fes_from_prob.py $f > $outfile
    head -1 $outfile |awk -v time=$t '{print time,$4}' >> fes_deltaF.${i}.data
  done
  cp `ls -v fes_running.${i}/* |tail -1` fes_running.${i}.data
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
