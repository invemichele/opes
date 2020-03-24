#! /usr/bin/env python3

### Create a reweighted FES ###

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#toggles
sigma=0.03
pace_to_time=2.5
print_stride=400

#setup
Kb=0.0083144621 #kj/mol
kbt=1
grid_min=-3
grid_max=3
grid_bin=100
transition_s=0
cv_grid=np.linspace(grid_min,grid_max,grid_bin)

#parser
parser = argparse.ArgumentParser(description='reweight')
parser.add_argument('-r',dest='replica',type=int,required=True,help='replica number')
parser.add_argument('-b',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-t',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('-f',dest='flip',action='store_true',default=False,required=False,help='flip time')
args = parser.parse_args()
wk='.'+str(args.replica)
print('  replica: '+wk)
bck=args.bck
if bck:
  print('  backup: '+bck)
tran=args.tran
flip=args.flip
sub_dir=''
if tran and flip:
  sys.exit(' choose either tran of flip')
if tran:
  sub_dir='tran'+str(tran)+'/'
  print('  tran=',tran)
if flip:
  sub_dir='tran-1/'
  print('  flip')
sub_dir='rct-'+sub_dir

#get colvar
filename=bck+'Colvar'+wk+'.data'
cv_col=1
bias_col=3
data=pd.read_csv(filename,sep='\s+',comment='#',header=None,usecols=[cv_col,bias_col])
filename=bck+'rct'+wk+'.data'
rct_col=1
data_rct=pd.read_csv(filename,sep='\s+',comment='#',header=None,usecols=[rct_col])
if flip:
  data=data.iloc[::-1]
  data_rct=data_rct.iloc[::-1]
cv=np.array(data.iloc[:,0])
bias=np.array(data.iloc[:,1])-np.array(data_rct.iloc[:,0])
del data
del data_rct

#output files
file_ext='.data'
fes_running_file='FES_rew'
head='cv_bin  fes'
current_fes_running=sub_dir+fes_running_file+wk+'/'+fes_running_file+'.t-%d'+file_ext
#create_dir='bck.meup.sh {0}; mkdir -p {0}'
create_dir='mkdir -p {0}'
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
  deltaF[n]=np.log((np.exp(-fes[cv_grid<transition_s]/kbt)).sum()/(np.exp(-fes[cv_grid>transition_s]/kbt)).sum())

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
if flip:
  head+=' # flip'
  time-=time[0]
  time=time[::-1]
#cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
#cmd.wait()
np.savetxt(filename,np.c_[time,deltaF],header=head,fmt='%14.9f')
if not tran:
  fes=-kbt*np.log(prob/max(prob))
  filename='FES_rew'+file_ext
  head='cv_bin  fes'
  #cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
  #cmd.wait()
  np.savetxt(filename,np.c_[cv_grid,fes],header=head,fmt='%14.9f')

