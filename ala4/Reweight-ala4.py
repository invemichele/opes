#! /usr/bin/env python3

### Create a reweighted FES ###

import sys
import numpy as np
import pandas as pd
import subprocess

#toggles
bck=''
#bck='bck.0.'
sigma=0.12 #taken from sigma used at the end of the OPES run
transient=0
if len(sys.argv)>1:
  transient=int(sys.argv[1])
  print('  using transient='+str(transient),file=sys.stderr)
#always skip comment rows
transient+=13

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
filename=bck+'Colvar.data'
x_col=2
y_col=3
bias_col=7
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[x_col,y_col,bias_col],skiprows=transient)
cv_x=np.array(data.ix[:,x_col])
cv_y=np.array(data.ix[:,y_col])
bias=np.array(data.ix[:,bias_col])
del data

#build fes
#basinA=0
#basinB=0
max_prob=0
prob=np.zeros((grid_bin,grid_bin))
for i in range(grid_bin):
  print('    working... {:.0%}'.format(i/grid_bin),end='\r',file=sys.stderr)
  for j in range(grid_bin):
    dx=np.absolute(x[i,j]-cv_x)
    dy=np.absolute(y[i,j]-cv_y)
    arg2=(np.minimum(dx,period-dx)/sigma)**2+(np.minimum(dy,period-dy)/sigma)**2
    prob[i,j]=np.sum(np.exp(bias/kbt)*np.exp(-0.5*arg2))
    if prob[i,j]>max_prob:
      max_prob=prob[i,j]
#    if x[i,j]>0 and x[i,j]<2.3:
#      basinB+=prob[i,j]
#    else:
#      basinA+=prob[i,j]

#print out
print('    printing...    ',end='\r',file=sys.stderr)
print('#cv_grid  fes')
#print('#cv_grid  fes #deltaF= ',kbt*np.log(basinA/basinB))
#print('#\n'*8,end='')
for i in range(grid_bin):
  for j in range(grid_bin):
    print(x[i,j],y[i,j],-kbt*np.log(prob[i,j]/max_prob))
  print('')

