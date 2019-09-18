#! /bin/bash

bck=''
#bck='bck.0.'
out=transition_time.data
for dir in `ls -d *-model2`
do
  cd $dir
  bck.meup.sh $out
  echo "#two_basins three_basins four_basins round_trip" > $out
  for f in `ls -v rep_*/${bck}Colvar.*`
  do 
    echo -n $f
    awk -v out=$out 'NR>1{
      if($3<-1.5){basin=1; if(a==0){a=1; printf " a";} if(end!=0 && end!=basin){printf "\n"; printf $1"\n" >> out; exit;} }
      if($3>1.5) {basin=2; if(c==0){c=1; printf " c";} if(end!=0 && end!=basin){printf "\n"; printf $1"\n" >> out; exit;} }
      if($2<-1.5){basin=3; if(b==0){b=1; printf " b";} if(end!=0 && end!=basin){printf "\n"; printf $1"\n" >> out; exit;} }
      if($2>1.5) {basin=4; if(d==0){d=1; printf " d";} if(end!=0 && end!=basin){printf "\n"; printf $1"\n" >> out; exit;} }
      new=a+b+c+d;
      if (new>old && new!=1) printf $1" " >> out; 
      if (end==0 && new==4) end=basin;
      old=new;
      }' $f
  done
  #./Median.py
  cd ..
done
