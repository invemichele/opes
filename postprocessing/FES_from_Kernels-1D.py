#! /usr/bin/env python3

### Get the running FES estimate used by OPES, 1D only ###
# similar to plumed sum_hills

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

# Manually change if necessary
recursive_merge=True

#parser
parser = argparse.ArgumentParser(description='get the running FES estimate used by OPES, 1D only')
parser.add_argument('--kernels',dest='filename',type=str,default='KERNELS',help='the kernels file name')
parser.add_argument('--kt',dest='kbt',type=float,required=True,help='the temperature in energy units')
parser.add_argument('--angle',dest='angle',action='store_true',default=False,help='the cv is an angle in the range [-pi,pi]')
parser.add_argument('--min',dest='grid_min',type=str,required=False,help='lower bound for the grid')
parser.add_argument('--max',dest='grid_max',type=str,required=False,help='upper bound for the grid')
parser.add_argument('--bin',dest='grid_bin',type=int,default=100,help='number of bins for the grid')
parser.add_argument('--stride',dest='stride',type=int,default=0,help='how often to print partial fes')
parser.add_argument('--mintozero',dest='mintozero',action='store_true',default=False,help='shift the minimum to zero')
parser.add_argument('--outfile',dest='outfile',type=str,default='fes.dat',help='grid spacing, alternative to number of bins')
args = parser.parse_args()
#parsing
filename=args.filename
kbt=args.kbt
mintozero=args.mintozero

#get kernels
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
for n in range(3):
  line=f.readline() #for safety skip some lines
line2=f.readline() #second line
pace_to_time=(float(line2.split()[0])-float(line.split()[0]))
f.close()
cutoff2=cutoff**2
val_at_cutoff=np.exp(-0.5*cutoff2)
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3])
center=np.array(data.iloc[:,0])
sigma=np.array(data.iloc[:,1])
height=np.array(data.iloc[:,2])
del data
print('  all data loaded')

#set grid
period=0
if args.angle:
  if (args.grid_min is not None or args.grid_max is not None):
    sys.exit('do not set min and max if variable is an angle')
  grid_min=-np.pi
  grid_max=np.pi
  period=2*np.pi
if args.grid_min is None:
  grid_min=min(center)
else:
  grid_min=args.grid_min
if args.grid_max is None:
  grid_max=max(center)
else:
  grid_max=args.grid_max
grid_bin=args.grid_bin
cv_grid=np.linspace(grid_min,grid_max,grid_bin)

#output files
head='cv  fes'
print_stride=args.stride
if print_stride==0:
  print_stride=len(center)+1
outfile=args.outfile
if print_stride<=len(center):
  file_ext=outfile.split('.')[-1]
  if len(file_ext)>1:
    file_ext='.'+file_ext
  current_fes_running='fes_running/'+outfile[:-len(file_ext)]+'.t-%d'+file_ext
  cmd=subprocess.Popen('mkdir fes_running',shell=True)
  cmd.wait()

# useful functions
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
    if period==0:
      arg=((cv_grid[x]-c)/s)**2
    else:
      dx=abs(cv_grid[x]-c)
      arg=(min(dx,period-dx)/s)**2
    if arg<cutoff2:
      delta[x]=h*(np.exp(-0.5*arg)-val_at_cutoff)
  return delta

def build_fes():
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
  norm=1
  if mintozero:
    norm=max(prob)
  return -kbt*np.log(prob/norm)

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
    z_fes=build_fes()
    np.savetxt(current_fes_running%((i+1)*pace_to_time),np.c_[cv_grid,z_fes],header=head,fmt='%14.9f')
print(' total kernels read from file: %d'%len(center))
print(' total kernels in compressed FES: %d'%len(z_center))

#build and print final
z_fes=build_fes()
np.savetxt(outfile,np.c_[cv_grid,z_fes],header=head,fmt='%14.9f')

