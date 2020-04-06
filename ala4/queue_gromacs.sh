#!/bin/bash

# Job Settings
jname=${PWD##*/}
ncore=1
max_t=24:00 #h:min
part=''
#part='-R "select[model==XeonGold_6150]"'
#singleton="-d singleton"
singleton=""

#to run locally
host=$HOSTNAME
[ $# -eq 1 ] && host=$1

# Commands
nsteps="-nsteps 10000000"
#nsteps=""
exe="`which gmx_mpi` mdrun"
Filename=alanine
Inputname=input
NumWalkers=`ls ${Inputname}*.tpr |wc -l`
ncore=$[ncore*NumWalkers]
WalkerPrefix="."
optWalkers="-multi $NumWalkers"
if [ $NumWalkers -eq 1 ]
then
  WalkerPrefix=""
  optWalkers=""
fi
Project=${Filename}${WalkerPrefix}
Run_file=${Inputname}${WalkerPrefix}.tpr
outfile=${Filename}.out
max_h=`python <<< "print('%g'%(${max_t%:*}+${max_t#*:}/60-0.05))"`

# Prepare Submission
bck.meup.sh -i $outfile
res=""
#cpt_files=`ls ${Project}*.cpt |wc -l`
#if [ $cpt_files -gt 0 ]
#then
#  res="-cpi ${Project}.cpt -append"
#  bck.meup.sh -i ${Project}*.gro > $outfile
#else
#  bck.meup.sh -i ${Project}* > $outfile
#fi
bck.meup.sh -i ${Project}* > $outfile

#mpi_cmd="$exe -maxh $max_h -s $Run_file -deffnm $Project $optWalkers $nsteps -ntomp 1"
mpi_cmd="$exe -plumed plumed.dat -maxh $max_h -s $Run_file -deffnm $Project $optWalkers $nsteps -ntomp 1 $res"
#extra_cmd="$0"
extra_cmd="../analyze-50.sh"

### if euler ###
if [ ${host:0:3} == "eu-" ]
then
  cmd="mpirun ${mpi_cmd}"
  if [ ! -z "$extra_cmd" ]
  then
    cmd="${cmd}; bsub -w \"done(${jname})\" -J after$jname -o $outfile $extra_cmd"
  fi
  submit="bsub -o $outfile -J $jname -n $ncore -W $max_t $part $cmd"
  echo -e " euler submission:\n$submit" |tee -a $outfile
### if daint ###
elif [ ${host:0:5} == "daint" ]
then
  hypt="--hint=nomultithread" #avoid hyperthreading
  sb=_sbatchme.sh
  echo -e "#!/bin/bash\nsrun $hypt ${mpi_cmd}\n${extra_cmd}" > $sb
  submit="sbatch -C mc -o $outfile -J $jname -n $ncore -t $max_t:00 $part $singleton $sb"
  echo -e " daint submission:\n$submit\n $sb:" |tee -a $outfile
  cat $sb |tee -a $outfile
### if workstation ###
else
  if [ $ncore -gt 8 ]
  then
    ncore=8
  fi
  submit="time mpirun -np $ncore ${mpi_cmd}"
  echo -e " workstation submission:\n$submit\n$extra_cmd" |tee -a $outfile
  eval "$submit &>> $outfile"
  submit="$extra_cmd" # &>> $outfile"
fi

# Actual Submission
eval $submit
