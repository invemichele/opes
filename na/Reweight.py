#! /usr/bin/env python3

### Create a reweighted FES ###

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#toggles
sigma=0.003
pace_to_time=1
print_stride=2000
rct_col=-1 #4

#setup
Kb=0.0083144621 #kj/mol
kbt=Kb*350
transition_s=1.99
grid_min=1.5
grid_max=2.7
grid_bin=100
cv_grid=np.linspace(grid_min,grid_max,grid_bin)

#parser
parser = argparse.ArgumentParser(description='reweight')
parser.add_argument('-r',dest='replica',type=int,default=-1,required=False,help='replica number')
parser.add_argument('-b',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-t',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
args = parser.parse_args()
wk=''
if args.replica!=-1:
  wk='.'+str(args.replica)
  print('  replica: '+wk)
bck=args.bck
if bck:
  print('  backup: '+bck)
tran=args.tran
if tran:
  sub_dir='tran'+str(tran)+'/'
  print('  tran=',tran)

#get colvar
filename=bck+'Colvar'+wk+'.data'
cv_col=1
bias_col=3
if rct_col>0:
  data=pd.read_csv(filename,sep='\s+',comment='#',header=None,usecols=[cv_col,bias_col,rct_col])
  bias=np.array(data.ix[:,bias_col])-np.array(data.ix[:,rct_col])
else:
  data=pd.read_csv(filename,sep='\s+',comment='#',header=None,usecols=[cv_col,bias_col])
  bias=np.array(data.ix[:,bias_col])
cv=np.array(data.ix[:,cv_col])
del data

#output files
file_ext='.data'
fes_running_file='FES_rew'
head='cv_bin  fes'
current_fes_running=sub_dir+fes_running_file+wk+'/'+fes_running_file+'.t-%d'+file_ext
create_dir='bck.meup.sh {0}; mkdir -p {0}'
cmd=subprocess.Popen(create_dir.format(sub_dir+fes_running_file+wk),shell=True)
cmd.wait()


n_tot=int(len(cv)/print_stride)
time=np.zeros(n_tot)
deltaF=np.zeros(n_tot)
prob=np.zeros(grid_bin)
def print_fes(it):
  fes=-kbt*np.log(prob/max(prob))
  np.savetxt(current_fes_running%((it+1)*pace_to_time),np.c_[cv_grid,fes],header=head,fmt='%14.9f')
  n=int(it/print_stride)
  time[n]=(it+1)*pace_to_time
  deltaF[n]=-np.log((np.exp(-fes[cv_grid<transition_s]/kbt)).sum()/(np.exp(-fes[cv_grid>transition_s]/kbt)).sum())

#build fes
for i in range(int(tran/pace_to_time),len(cv)):
  print('    working... {:.0%}'.format(i/len(cv)),end='\r')
  arg=(cv_grid-cv[i])/sigma
  prob+=np.exp(bias[i]/kbt)*np.exp(-0.5*arg**2)
  if (i+1)%print_stride==0:
    print_fes(i)

#output files
file_ext=wk+'.data'
filename=sub_dir+'fes_deltaF.rew'+file_ext
head='time  deltaF'
cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
cmd.wait()
np.savetxt(filename,np.c_[time,deltaF],header=head,fmt='%14.9f')
if not tran:
  fes=-kbt*np.log(prob/max(prob))
  filename='FES_rew'+file_ext
  head='cv_bin  fes'
  cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
  cmd.wait()
  np.savetxt(filename,np.c_[cv_grid,fes],header=head,fmt='%14.9f')

