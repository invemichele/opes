#! /usr/bin/env python3

### Get the running fes estimate used by OPES ###
# see ../../postprocessing/FES_from_Kernels-1D.py for a more updated version

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#toggles
print_stride=400
recursive_merge=True
print_compressed_kernels=False

#setup
Kb=0.0083144621 #kj/mol
kbt=1
transition_s=0
grid_min=-3
grid_max=3
grid_bin=100
cv_grid=np.linspace(grid_min,grid_max,grid_bin)

#parser
parser = argparse.ArgumentParser(description='get the FES running estimate')
parser.add_argument('-r',dest='replica',type=int,required=True,help='replica number')
parser.add_argument('-b',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
args = parser.parse_args()
wk='.'+str(args.replica)
print('  replica: '+wk)
bck=args.bck
if bck:
  print('  backup: '+bck)

##get walkers num
#cmd=subprocess.Popen('ls Colvar.* |wc -l',shell=True,stdout=subprocess.PIPE)
#output=cmd.communicate()
#n_walkers=int(output[0])
#if(n_walkers>0):
#  print('  - n_walkers found: %d'%n_walkers)
#  print_stride*=n_walkers
#else:
#  n_walkers=1
n_walkers=1

#get kernels
filename=bck+'Kernels'+wk+'.data'
f=open(filename)
line=f.readline() #header
if len(line.split())!=7:
  sys.exit('  something is wrong with file '+filename)
line=f.readline() #biasfactor
line=f.readline() #epsilon
epsilon=float(line.split()[-1])
line=f.readline() #cutoff
cutoff=float(line.split()[-1])
line=f.readline() #threshold
threshold=float(line.split()[-1])
for n in range(n_walkers):
  line=f.readline() #first line
line2=f.readline() #second line
pace_to_time=(float(line2.split()[0])-float(line.split()[0]))/n_walkers
f.close()
cutoff2=cutoff**2
val_at_cutoff=np.exp(-0.5*cutoff2)
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3])
center=np.array(data.iloc[:,0])
sigma=np.array(data.iloc[:,1])
height=np.array(data.iloc[:,2])
del data
print('  all data loaded',end='\r')

#output files
file_ext='.data'
fes_running_file='fes_running'
head='cv_bin  fes'
current_fes_running=fes_running_file+wk+'/'+fes_running_file+'.t-%d'+file_ext
#create_dir='bck.meup.sh {0}; mkdir -p {0}'
create_dir='mkdir -p {0}'
cmd=subprocess.Popen(create_dir.format(fes_running_file+wk),shell=True)
cmd.wait()

# useful functions
n_tot=int(len(center)/print_stride)
time=np.zeros(n_tot)
deltaF=np.zeros(n_tot)

z_center=[center[0]]
z_sigma=[sigma[0]]
z_height=[height[0]]

def get_merge_candidate(c,self):
  min_dist=threshold
  min_j=-1
  for j in range(len(z_center)):
    if j==self:
      continue
    dist=abs(c-z_center[j])/z_sigma[j]
    if dist<min_dist:
      min_j=j
  return min_j

def merge(j,m_height,m_center,m_sigma):
    h=z_height[j]+m_height
    c=(z_height[j]*z_center[j]+m_height*m_center)/h
    s2=(z_height[j]*(z_sigma[j]**2+z_center[j]**2)+m_height*(m_sigma**2+m_center**2))/h-c**2
    z_height[j]=h
    z_center[j]=c
    z_sigma[j]=np.sqrt(s2)

def delta_kernel(c,s,h):
  delta=np.zeros(len(cv_grid))
  for x in range(len(cv_grid)):
    arg=((cv_grid[x]-c)/s)**2
    if arg<cutoff2:
      delta[x]=h*(np.exp(-0.5*arg)-val_at_cutoff)
  return delta

def print_fes(it):
  prob=np.zeros(grid_bin)
  for j in range(len(z_center)):
    prob+=delta_kernel(z_center[j],z_sigma[j],z_height[j])
  Zed=0
  for j in range(len(z_center)):
    for jj in range(len(z_center)):
      arg=((z_center[jj]-z_center[j])/z_sigma[j])**2
      if arg<cutoff2:
        Zed+=z_height[j]*(np.exp(-0.5*arg)-val_at_cutoff)
  Zed/=len(z_center)
  prob=prob/Zed+epsilon
  z_fes=-kbt*np.log(prob/max(prob))
  np.savetxt(current_fes_running%((it+1)*pace_to_time),np.c_[cv_grid,z_fes],header=head,fmt='%14.9f')
  #calc other stuff
  n=int(it/print_stride)
  time[n]=(it+1)*pace_to_time
  deltaF[n]=np.log((np.exp(-z_fes[cv_grid<transition_s]/kbt)).sum()/(np.exp(-z_fes[cv_grid>transition_s]/kbt)).sum())

# compression
for i in range(1,len(center)):
  print('    working... {:.0%}'.format(i/len(center)),end='\r')
  j=get_merge_candidate(center[i],-1)
  if j>=0:
    merge(j,height[i],center[i],sigma[i])
    if recursive_merge:
      jj=get_merge_candidate(z_center[j],j)
      while jj>=0:
        merge(jj,z_height[j],z_center[j],z_sigma[j])
        z_height.pop(j)
        z_center.pop(j)
        z_sigma.pop(j)
        if j<jj:
          j=jj-1
        else:
          j=jj
        jj=get_merge_candidate(z_center[j],j)
  else:
    z_center.append(center[i])
    z_sigma.append(sigma[i])
    z_height.append(height[i])
  if (i+1)%print_stride==0:
    print_fes(i)
print(' total kernels in compressed fes: %d'%len(z_center))
print(' compression rate: %g'% ((len(center)-len(z_center))/len(center)*100) )

#print deltaF
filename='fes_deltaF'+wk+file_ext
#backup='bck.meup.sh '+filename
#cmd=subprocess.Popen(backup,shell=True)
#cmd.wait()
head='time  deltaF #transition_s=%g'%transition_s
np.savetxt(filename,np.c_[time,deltaF],header=head,fmt='%14.9f')

#build and print final
#backup='bck.meup.sh '+fes_running_file+wk+file_ext
#cmd=subprocess.Popen(backup,shell=True)
#cmd.wait()
prob=np.zeros(grid_bin)
for i in range(len(z_center)):
  prob+=delta_kernel(z_center[i],z_sigma[i],z_height[i])
Zed=0
for j in range(len(z_center)):
  for jj in range(len(z_center)):
    arg=((z_center[jj]-z_center[j])/z_sigma[j])**2
    if arg<cutoff2:
      Zed+=z_height[j]*(np.exp(-0.5*arg)-val_at_cutoff)
Zed/=len(z_center)
prob=prob/Zed+epsilon
z_fes=-kbt*np.log(prob/max(prob))
head='cv  z_fes #threshold=%g'%threshold
np.savetxt(fes_running_file+wk+file_ext,np.c_[cv_grid,z_fes],header=head,fmt='%14.9f')

#print final kernels
if print_compressed_kernels:
  filename='compressed-Kernels'+wk+'.data'
  #backup='bck.meup.sh -i '+filename
  #cmd=subprocess.Popen(backup,shell=True)
  #cmd.wait()
  head='center  sigma  height #threshold=%g'%threshold
  np.savetxt(filename,np.c_[z_center,z_sigma,z_height],header=head,fmt='%14.9f')
