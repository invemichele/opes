#! /bin/bash

bck=''
#bck='bck.0.'

for dir in `ls -d rep_*`
do
  cd $dir
  n_rep=`ls Colvar.* |wc -l`
  echo " number of replicas: $n_rep"

  if [ $n_rep -lt 2 ]
  then
    echo " replicas not found"
    exit
  fi

  for i in `seq 0 $[n_rep-1]`
  do
    echo -e "\n\n $dir/$i \n"
    ../Build-fes_running.py -r $i
    ../Reweight-multi.py -r $i -f
  done
  cd ..
done

make_stats() {
  bck.meup.sh -i $2
  paste $1 |\
    awk 'BEGIN{print "#average std_dev"}
         NR>1{
           av=0; av2=0; 
           for (i=2; i<=NF; i+=2) {av+=$i; av2+=$i^2;};
           av*=2./NF; av2*=2./NF; 
           print $1,av,sqrt(av2-av^2)
         }' > $2
}

make_stats "rep_?/fes_deltaF.?.data" stats-fes_deltaF.data
make_stats "rep_?/tran-1/fes_deltaF.rew.?.data" flip-stats-fes_deltaF.rew.data
