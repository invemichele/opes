#!/bin/bash

# Job Settings
jname=m-${PWD##*/}
ncore=10
max_t=4:00 #h:min
#part=parrinello_compute #for monch only
part=express_compute
#logfile=log.md
outfile=log.out

#to run locally
host=$HOSTNAME
[ $# -eq 1 ] && host=$1

# Commands
mpi_cmd="plumed ves_md_linearexpansion md_input.dat"
extra_cmd=""

# Prepare Submission
bck.meup.sh -i $outfile
#bck.meup.sh -v $logfile |& tee -a $outfile
### if euler ###
if [ ${host:0:3} == "eu-" ]
then
  cmd="mpirun ${mpi_cmd}; ${extra_cmd}"
  submit="bsub -o $outfile -J $jname -n $ncore -W $max_t $cmd"
  echo -e " euler submission:\n$submit" |tee -a $outfile
### if monch ###
elif [ ${host:0:5} == "monch" ]
then
  hypt="--ntasks-per-core 1" #set to 1 for avoiding hyperthreading
  sb=_sbatchme.sh
  echo -e "#!/bin/bash\nmpirun -np $ncore ${mpi_cmd}\n${extra_cmd}" > $sb
  submit="sbatch -o $outfile -J $jname -n $ncore $hypt -t $max_t:00 -p $part $sb"
  echo -e " monch submission:\n$submit\n $sb:" |tee -a $outfile
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
$submit
#eval $submit
