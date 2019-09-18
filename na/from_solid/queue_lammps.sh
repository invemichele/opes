#!/bin/bash

# Job Settings
jname=${PWD##*/}
max_t=120:00 #h:min
ncore=10
nrep=4
part=''
#part='-R "select[model==XeonGold_6150]"'

#to run locally
host=$HOSTNAME
[ $# -eq 1 ] && host=$1

# Commands
logfile=log.lammps
outfile=log.out
exe=`which lmp_mpi`
max_sec=$((${max_t%:*}*3600+${max_t#*:}*60-120)) #two min less
#out_op="-echo screen -screen $logfile -log none"
out_op="-screen none"
mw=''
if [ $nrep -gt 1 ]
then
  mw="-p ${nrep}x${ncore}"
fi
ncore=$[nrep*ncore]

#mpi_cmd="$exe $mw -var max_time $max_sec -in input-restart.lmp $out_op"
mpi_cmd="$exe $mw -var max_time $max_sec -in input.lmp $out_op"
extra_cmd="thermo_extractor.sh"

# Prepare Submission
bck.meup.sh -i $outfile
bck.meup.sh -v ${logfile}* |& tee -a $outfile
### if euler ###
if [ ${HOSTNAME:0:3} == "eu-" ]
then
  cmd="mpirun ${mpi_cmd}; ${extra_cmd}"
  submit="bsub -o $outfile -J $jname -n $ncore -W $max_t $part $cmd"
  echo -e " euler submission:\n$submit" |tee -a $outfile
### if monch ###
elif [ ${HOSTNAME:0:5} == "monch" ]
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
