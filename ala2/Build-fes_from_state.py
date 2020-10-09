#! /usr/bin/env python3

### Create a reweighted FES from a state file ###
# see ../postprocessing/FES_from_State-2D.py for a more updated version

import sys
import numpy as np
import pandas as pd

if len(sys.argv)<2:
  sys.exit(' I need the filename of the compressed kernels to build the fes')

#setup
Kb=0.0083144621 #kj/mol
kbt=Kb*300
grid_min=-np.pi
grid_max=np.pi
grid_bin=100
cv_grid=np.linspace(grid_min,grid_max,grid_bin)
x,y=np.meshgrid(cv_grid,cv_grid)
period=grid_max-grid_min

#get kernels
filename=sys.argv[1]
f=open(filename)
line=f.readline()
if len(line.split())!=8:
  sys.exit('  something is wrong with file'+filename)
line=f.readline() #action
line=f.readline() #biasfactor
line=f.readline() #epsilon
epsilon=float(line.split()[-1])
line=f.readline() #kernel_cutoff
cutoff=float(line.split()[-1])
val_at_cutoff=np.exp(-0.5*cutoff**2)
line=f.readline() #compression_threshold
line=f.readline() #zed
Zed=float(line.split()[-1])
f.close()
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3,4,5])
center_x=np.array(data.iloc[:,0])
center_y=np.array(data.iloc[:,1])
sigma_x=np.array(data.iloc[:,2])
sigma_y=np.array(data.iloc[:,3])
height=np.array(data.iloc[:,4])
del data

#calculate
basinA=0
basinB=0
max_prob=0
prob=np.zeros((grid_bin,grid_bin))
for i in range(grid_bin):
  print('    working... {:.0%}'.format(i/grid_bin),end='\r',file=sys.stderr)
  for j in range(grid_bin):
    dx=np.absolute(x[i,j]-center_x)
    dy=np.absolute(y[i,j]-center_y)
    arg2=(np.minimum(dx,period-dx)/sigma_x)**2+(np.minimum(dy,period-dy)/sigma_y)**2
    prob[i,j]=np.sum(height*(np.maximum(np.exp(-0.5*arg2)-val_at_cutoff,0)))
    prob[i,j]=prob[i,j]/Zed+epsilon
    if prob[i,j]>max_prob:
      max_prob=prob[i,j]
    if x[i,j]>0 and x[i,j]<2.3:
      basinB+=prob[i,j]
    else:
      basinA+=prob[i,j]

#print out
print('    printing...    ',end='\r',file=sys.stderr)
print('#cv_grid  fes #deltaF= ',kbt*np.log(basinA/basinB))
for i in range(grid_bin):
  for j in range(grid_bin):
    print(x[i,j],y[i,j],-kbt*np.log(prob[i,j]/max_prob))
  print('')
