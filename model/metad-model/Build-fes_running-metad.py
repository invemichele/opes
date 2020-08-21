#! /usr/bin/env python3

### home-made sum_hills + real rct ###

import sys
import numpy as np
import pandas as pd
import linecache
import subprocess
import argparse

#toggles
pace_to_time=2.5
print_stride=400
DP2CUTOFF=6.25 #same used by METAD
interval=2.9 #needed if using INTERVAL keyword

#setup
Kb=0.0083144621 #kj/mol
kbt=1
grid_min=-3
grid_max=3
grid_bin=100
transition_s=0
cv_grid=np.linspace(grid_min,grid_max,grid_bin)
cv_grid_interval=np.array(cv_grid)
for x in range(len(cv_grid_interval)):
  if cv_grid_interval[x]>interval:
    cv_grid_interval[x]=interval
  elif cv_grid_interval[x]<-interval:
    cv_grid_interval[x]=-interval

#parser
parser = argparse.ArgumentParser(description='reweight')
parser.add_argument('-r',dest='replica',type=int,required=True,help='replica number')
parser.add_argument('-b',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
args = parser.parse_args()
wk='.'+str(args.replica)
print('  replica: '+wk)
bck=args.bck
if bck:
  print('  backup: '+bck)

#get hills
filename=bck+'Hills'+wk+'.data'
one_hills_line=linecache.getline(filename,100)
if len(one_hills_line.split()) != 5:
  sys.exit(' ERROR: hills file not compatible')
b_sigma=float(one_hills_line.split()[2])
gamma=float(one_hills_line.split()[4])
inv_gamma=0 #non well-tempered case
if gamma!=-1:
  inv_gamma=1/gamma
data=pd.read_csv(filename,sep='\s+',comment='#',header=None,usecols=[0,1,3])
b_time=np.array(data.iloc[:,0])
b_center=np.array(data.iloc[:,1])
b_height=np.array(data.iloc[:,2])
del data

#get true fes
data=pd.read_csv('../../ref_fes-model.data',sep='\s+',comment='#',header=None,usecols=[1])
true_fes=np.array(data.iloc[:,0])
del data
Z0=(np.exp(-true_fes/kbt)).sum()
true_fes+=kbt*np.log(Z0) #normalization

#output files
file_ext='.data'
fes_running_file='fes_running'
head='cv_bin  fes'
current_fes_running=fes_running_file+wk+'/'+fes_running_file+'.t-%d'+file_ext
#create_dir='bck.meup.sh {0}; mkdir -p {0}'
create_dir='mkdir -p {0}'
cmd=subprocess.Popen(create_dir.format(fes_running_file+wk),shell=True)
cmd.wait()

#variables
n_tot=int(len(b_center)/print_stride)
time=np.zeros(n_tot)
deltaF=np.zeros(n_tot)
bias=np.zeros(grid_bin)
rct=np.zeros(len(b_center))

def print_fes(it):
  fes=-1/(1-inv_gamma)*bias
  fes-=min(fes)
  np.savetxt(current_fes_running%((it+1)*pace_to_time),np.c_[cv_grid,fes],header=head,fmt='%14.9f')
  n=int(it/print_stride)
  time[n]=(it+1)*pace_to_time
  deltaF[n]=np.log((np.exp(-fes[cv_grid<transition_s]/kbt)).sum()/(np.exp(-fes[cv_grid>transition_s]/kbt)).sum())

#build fes
for i in range(len(b_center)):
  print('    working... {:.0%}'.format(i/len(b_center)),end='\r')
  dp2=0.5*((cv_grid_interval-b_center[i])/b_sigma)**2
  for x in range(len(dp2)):
    if dp2[x]>DP2CUTOFF:
      dp2[x]=np.inf
  bias+=(1-inv_gamma)*b_height[i]*np.exp(-dp2)
  rct[i]=-kbt*np.log((np.exp(-(true_fes+bias)/kbt)).sum())
  if (i+1)%print_stride==0:
    print_fes(i)

#output files
file_ext=wk+'.data'
filename='fes_deltaF'+file_ext
head='time  deltaF'
#cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
#cmd.wait()
np.savetxt(filename,np.c_[time,deltaF],header=head,fmt='%14.9f')

filename='rct'+file_ext
head='# time  rct\n0 0'
#cmd=subprocess.Popen('bck.meup.sh -i '+filename,shell=True)
#cmd.wait()
np.savetxt(filename,np.c_[b_time,rct],header=head,comments='',fmt='%14.9f')
