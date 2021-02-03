#! /usr/bin/env python3

### Generate an OPES STATE file from a KERNELS file ###
# For postprocessing only, do not use for restarting a simulation
# (the idea is to fake a restart with the plumed driver and dump the OPES state)

# WARNING: requires small modifications that will be added to PLUMED: https://github.com/invemichele/plumed2/commit/7be6829ce33ab22ce3ec0ff7ef1788ec638eef75

import sys
import argparse
import subprocess
plumed_exe='plumed'

#parser
parser = argparse.ArgumentParser(description='Generate an OPES STATE file from a KERNELS file, so that it can be used with FES_from_State.py. DO NOT use the obtained STATE file for restart')
parser.add_argument('--kernels','-f',dest='filename',type=str,default='KERNELS',help='the kernels file name, with the deposited kernels')
parser.add_argument('--outfile','-o',dest='outfile',type=str,default='STATE',help='name of the output file')
parser.add_argument('--keep_tmp',dest='keep_tmp',action='store_true',default=False,help='keep the temporary plumed file')
args = parser.parse_args()
#parsing
filename=args.filename
outfile=args.outfile

#get info
f=open(filename,'r')
line=f.readline() #fields
if line.split()[1]!='FIELDS':
  sys.exit(' no FIELDS found in file "'+filename+'"')
if len(line.split())<7:
  sys.exit(' not enough FIELDS found in file "'+filename+'"')
if (len(line.split())-5)%2!=0:
  sys.exit(' wrong number of FIELDS found in file "'+filename+'"')
ncv=int((len(line.split())-5)/2)
cvname=[]
for i in range(ncv):
  cvname.append(line.split()[3+i])
  if cvname[i].find('.')!=-1:
    sys.exit(' %s: unfourtunately, you must modify the KERNELS file and remove any "." from CVs names'%cvname[i])
line=f.readline() #action
if line.split()[3]=='OPES_METAD_kernels':
  action='OPES_METAD'
elif line.split()[3]=='OPES_METAD_EXPLORE_kernels':
  action='OPES_METAD_EXPLORE'
else:
  sys.exit(' this script only works with OPES_METAD or OPES_METAD_EXPLORE KERNELS files')
line=f.readline() #biasfactor
if line.split()[2]!='biasfactor':
  sys.exit(' biasfactor not found!')
biasfactor=line.split()[3]
line=f.readline() #epsilon
if line.split()[2]!='epsilon':
  sys.exit(' epsilon not found!')
epsilon=line.split()[3]
line=f.readline() #kernel_cutoff
if line.split()[2]!='kernel_cutoff':
  sys.exit(' kernel_cutoff not found!')
kernel_cutoff=line.split()[3]
line=f.readline() #compression_threshold
if line.split()[2]!='compression_threshold':
  sys.exit(' compression_threshold not found!')
compression_threshold=line.split()[3]
periodic=['NO']*ncv
line=f.readline()
while line.split()[0]=='#!':
  for i in range(ncv):
    if line.split()[2]=='min_'+cvname[i]:
      periodic[i]=line.split()[3]+','
      line=f.readline()
      if line.split()[2]!='max_'+cvname[i]:
        sys.exit(' periodic CVs should have both min and max value!')
      periodic[i]+=line.split()[3]
  line=f.readline()
f.close()

#create temporary plumed file
plumed_input='# vim:ft=plumed\n'
plumed_input+='RESTART\n'
plumed_input+='f: FIXEDATOM AT=0,0,0\n' #fake atom
plumed_input+='d: DISTANCE ATOMS=f,f\n' #unfourtunately the FAKE colvar has issues with PERIODIC
plumed_input+='COMMITTOR ARG=d BASIN_LL1=-1 BASIN_UL1=1\n' #this will kill the driver
for i in range(ncv):
  plumed_input+=cvname[i]+': COMBINE ARG=d PERIODIC='+periodic[i]+'\n' #recreate CVs label doesn't work if they have components!
plumed_input+='opes: '+action
plumed_input+=' TEMP=1 PACE=1 BARRIER=0'
plumed_input+=' ARG='+cvname[0]
for ii in range(1,ncv):
  plumed_input+=','+cvname[ii]
plumed_input+=' FILE='+filename
plumed_input+=' STATE_WFILE='+outfile
plumed_input+=' BIASFACTOR='+biasfactor
plumed_input+=' EPSILON='+epsilon
plumed_input+=' KERNEL_CUTOFF='+kernel_cutoff
plumed_input+=' COMPRESSION_THRESHOLD='+compression_threshold

tmp_plumed_file='tmp-plumed_driver.dat'
f=open(tmp_plumed_file,'w')
f.write(plumed_input)
f.close()

#run driver
cmd_string=plumed_exe+' driver --noatoms --plumed '+tmp_plumed_file
if not args.keep_tmp:
  cmd_string+='; rm '+tmp_plumed_file
cmd=subprocess.Popen(cmd_string,shell=True)
cmd.wait()

