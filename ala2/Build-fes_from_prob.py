#! /usr/bin/env python3

### Create a reweighted FES ###

import sys
import numpy as np
import pandas as pd
import subprocess

if len(sys.argv)<2:
  sys.exit(' I need the filename of the kernels to build the fes')

mode='EXPLORE'

if mode=='EXPLORE':
  biasfactor=10
else:
  biasfactor=1

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
line=f.readline() #cutoff
cutoff=float(line.split()[-1])
val_at_cutoff=np.exp(-0.5*cutoff**2)
#print('  cutoff=%g'%cutoff,file=sys.stderr)
#print('  val_at_cutoff=%g'%val_at_cutoff,file=sys.stderr)
line=f.readline() #epsilon
epsilon=float(line.split()[-1])
#print('  epsilon=%g'%epsilon,file=sys.stderr)
#line=f.readline() #norm
#norm=floar(line.split()[-1])
#print('  norm=%g'%norm,file=sys.stderr)
f.close()
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3,4,5])
center_x=np.array(data.ix[:,1])
center_y=np.array(data.ix[:,2])
sigma_x=np.array(data.ix[:,3])
sigma_y=np.array(data.ix[:,4])
height=np.array(data.ix[:,5])
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
    prob[i,j]=np.sum(height*np.maximum(np.exp(-0.5*arg2)-val_at_cutoff,0))
    prob[i,j]+=epsilon
    prob[i,j]=np.power(prob[i,j],biasfactor)
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
