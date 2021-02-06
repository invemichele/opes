#! /usr/bin/env python3

### Get the FES estimate from reweighting. 1D or 2D only ###
# uses a weighted kernel density estimation, so it requires the bandwidth sigma
# usage is similar to plumed sum_hills

import sys
import argparse
import numpy as np
import pandas as pd #much faster reading from file
#use_bck=False #requires the bck.meup.sh script
use_bck=True
if use_bck:
  import subprocess

print('')
### Parser stuff ###
parser = argparse.ArgumentParser(description='calculate the free energy surfase (FES) along the chosen collective variables (1 or 2) using a reweighted kernel density estimate')
# files
parser.add_argument('--colvar','-f',dest='filename',type=str,default='COLVAR',help='the COLVAR file name, with the collective variables and the bias')
parser.add_argument('--outfile','-o',dest='outfile',type=str,default='fes_rew.dat',help='name of the output file')
# compulsory
parser.add_argument('--sigma','-s',dest='sigma',type=str,required=True,help='the bandwidth for the kernel density estimation. Use e.g. the last value of sigma from an OPES_METAD simulation')
kbt_group = parser.add_mutually_exclusive_group(required=True)
kbt_group.add_argument('--kt',dest='kbt',type=float,help='the temperature in energy units')
kbt_group.add_argument('--temp',dest='temp',type=float,help='the temperature. Energy units is Kj/mol')
# input columns
parser.add_argument('--cv',dest='cv',type=str,default='2',help='the CVs to be used. Either by name or by column number, starting from 1')
parser.add_argument('--bias',dest='bias',type=str,default='.bias',help='the bias to be used. Either by name or by column number, starting from 1. Set to NO for nonweighted KDE')
# grid related
parser.add_argument('--min',dest='grid_min',type=str,required=False,help='lower bounds for the grid')
parser.add_argument('--max',dest='grid_max',type=str,required=False,help='upper bounds for the grid')
parser.add_argument('--bin',dest='grid_bin',type=str,default="100,100",help='number of bins for the grid')
# other options
parser.add_argument('--fmt',dest='fmt',type=str,default='% 12.6f',help='specify the output format')
parser.add_argument('--deltaFat',dest='deltaFat',type=float,required=False,help='calculate the free energy difference between left and right of given cv1 value')
parser.add_argument('--mintozero',dest='mintozero',action='store_true',default=False,help='shift the minimum to zero')
parser.add_argument('--der',dest='der',action='store_true',default=False,help='calculate also FES derivatives')
# some easy parsing
args=parser.parse_args()
if args.kbt is not None:
  kbt=args.kbt
else:
  kbt=args.temp*0.0083144621
fmt=args.fmt
calc_der=args.der
calc_deltaF=False
if args.deltaFat is not None:
  calc_deltaF=True
ts=args.deltaFat

### Get data ###
# get dim
dim=len(args.cv.split(','))
if dim==1:
  dim2=False
elif dim==2:
  dim2=True
else:
  sys.exit(' only 1D and 2D are supported')
# get cvs
f=open(args.filename,'r')
fields=f.readline().split()
if fields[1]!='FIELDS':
  sys.exit(' no FIELDS found in "%s"'%args.filename)
try:
  col_x=int(args.cv.split(',')[0])-1
  name_cv_x=fields[col_x+2]
except ValueError:
  col_x=-1
  name_cv_x=args.cv.split(',')[0]
  for i in range(len(fields)):
    if fields[i]==name_cv_x:
      col_x=i-2
  if col_x==-1:
    sys.exit(' cv "%s" not found'%name_cv_x)
  pass
if dim2:
  try:
    col_y=int(args.cv.split(',')[1])-1
    name_cv_y=fields[col_y+2]
  except ValueError:
    col_y=-1
    name_cv_y=args.cv.split(',')[1]
    for i in range(len(fields)):
      if fields[i]==name_cv_y:
        col_y=i-2
    if col_y==-1:
      sys.exit(' cv "%s" not found'%name_cv_y)
    pass
# get bias
if args.bias=='NO':
  col_bias=[]
else:
  try:
    col_bias=[int(col)-1 for col in args.bias.split(',')]
  except ValueError:
    col_bias=[]
    if args.bias=='.bias':
      for i in range(len(fields)):
        if fields[i].find('.bias')!=-1:
          col_bias.append(i-2)
    else:
      for j in range(len(args.bias.split(','))):
        for i in range(len(fields)):
          if fields[i]==args.bias.split(',')[j]:
            col_bias.append(i-2)
      if len(col_bias)!=len(args.bias.split(',')):
        sys.exit(' found %d matching bias, but %d were requested. Use columns number to avoid ambiguity'%(len(col_bias),len(args.bias.split(','))))
    pass
# get periodicity
period_x=0
period_y=0
line=f.readline().split()
while line[0]=='#!':
  if line[2]=='min_'+name_cv_x:
    if line[3]=='-pi':
      grid_min_x=-np.pi
    else:
      grid_min_x=float(line[3])
    line=f.readline().split()
    if line[2]!='max_'+name_cv_x:
      sys.exit(' min_%s was found, but not max_%s !'%(name_cv_x,name_cv_x))
    if line[3]=='pi':
      grid_max_x=np.pi
    else:
      grid_max_x=float(line[3])
    period_x=grid_max_x-grid_min_x
    if calc_der:
      sys.exit(' derivatives not supported with periodic CVs, remove --der option')
  if dim2 and line[2]=='min_'+name_cv_y:
    if line[3]=='-pi':
      grid_min_y=-np.pi
    else:
      grid_min_y=float(line[3])
    line=f.readline().split()
    if line[2]!='max_'+name_cv_y:
      sys.exit(' min_%s was found, but not max_%s !'%(name_cv_y,name_cv_y))
    if line[3]=='pi':
      grid_max_y=np.pi
    else:
      grid_max_y=float(line[3])
    period_y=grid_max_y-grid_min_y
    if calc_der:
      sys.exit(' derivatives not supported with periodic CVs, remove --der option')
  line=f.readline().split()
f.close()
# get sigma
sigma_x=float(args.sigma.split(',')[0])
if dim2:
  if len(args.sigma.split(','))!=2:
    sys.exit(' two comma-separated floats expected after --sigma')
  sigma_y=float(args.sigma.split(',')[1])
# read file
all_cols=[col_x]+col_bias
if dim2:
  all_cols=[col_x,col_y]+col_bias
data=pd.read_table(args.filename,dtype=float,sep='\s+',comment='#',header=None,usecols=all_cols)
cv_x=np.array(data.iloc[:,0])
it=1
if dim2:
  cv_y=np.array(data.iloc[:,1])
  it=2
if len(col_bias)==0:
  bias=np.zeros(len(cv_x))
else:
  bias=np.array(data.iloc[:,it])
  for i in range(1,len(col_bias)):
    it+=1
    bias+=np.array(data.iloc[:,it])
  bias/=kbt #dimensionless bias

### Prepare the grid ###
grid_bin_x=int(args.grid_bin.split(',')[0])
if period_x==0:
  grid_bin_x+=1 #same as plumed sum_hills
if args.grid_min is None:
  if period_x==0: #otherwise is already set
    grid_min_x=min(center_x)
else:
  if args.grid_min.split(',')[0]=='-pi':
    grid_min_x=-np.pi
  else:
    grid_min_x=float(args.grid_min.split(',')[0])
if args.grid_max is None:
  if period_x==0: #otherwise is already set
    grid_max_x=max(center_x)
else:
  if args.grid_max.split(',')[0]=='pi':
    grid_max_x=np.pi
  else:
    grid_max_x=float(args.grid_max.split(',')[0])
cv_grid_x=np.linspace(grid_min_x,grid_max_x,grid_bin_x)
if dim2:
  if len(args.grid_bin.split(','))!=2:
    sys.exit('two comma separated integers expected after --bin')
  grid_bin_y=int(args.grid_bin.split(',')[1])
  if period_y==0:
    grid_bin_y+=1 #same as plumed sum_hills
  if args.grid_min is None:
    if period_y==0: #otherwise is already set
      grid_min_y=min(center_y)
  else:
    if len(args.grid_min.split(','))!=2:
      sys.exit('two comma separated floats expected after --min')
    if args.grid_min.split(',')[1]=='-pi':
      grid_min_y=-np.pi
    else:
      grid_min_y=float(args.grid_min.split(',')[1])
  if args.grid_max is None:
    if period_y==0: #otherwise is already set
      grid_max_y=max(center_y)
  else:
    if len(args.grid_max.split(','))!=2:
      sys.exit('two comma separated floats expected after --max')
    if args.grid_max.split(',')[1]=='pi':
      grid_max_y=np.pi
    else:
      grid_max_y=float(args.grid_max.split(',')[1])
  cv_grid_y=np.linspace(grid_min_y,grid_max_y,grid_bin_y)
  x,y=np.meshgrid(cv_grid_x,cv_grid_y)
if calc_deltaF and (ts<=grid_min_x or ts>=grid_max_x):
  print(' +++ WARNING: the provided --deltaFat is out of the CV grid +++')
  calc_deltaF=False

### Calculate FES ###
# on single grid point
def calcFESpoint(point_x,point_y=None):
  if period_x==0:
    dist_x=(point_x-cv_x)/sigma_x
  else:
    dx=np.absolute(point_x-cv_x)
    dist_x=np.minimum(dx,period_x-dx)/sigma_x
  arg=bias-0.5*dist_x*dist_x
  if point_y is not None:
    if period_y==0:
      dist_y=(point_y-cv_y)/sigma_y
    else:
      dy=np.absolute(point_y-cv_y)
      dist_y=np.minimum(dy,period_y-dy)/sigma_y
    arg-=0.5*dist_y*dist_y
  if calc_der:
    arg_max=np.amax(arg)
    safe_kernels=np.exp(arg-arg_max)
    safe_prob=np.sum(safe_kernels)
    _fes=-kbt*(arg_max+np.log(safe_prob))
    _der_fes_x=-kbt*(np.sum(-dist_x/sigma_x*safe_kernels)/safe_prob)
    if point_y is None:
      return _fes,_der_fes_x
    else:
      _der_fes_y=-kbt*(np.sum(-dist_y/sigma_y*safe_kernels)/safe_prob)
      return _fes,_der_fes_x,_der_fes_y
  else:
    return -kbt*np.logaddexp.reduce(arg)
# loop over whole grid
if not dim2:
  fes=np.zeros(grid_bin_x)
  if calc_der:
    der_fes_x=np.zeros(grid_bin_x)
    for i in range(grid_bin_x):
      print('   working...  {:.0%} of '.format(i/grid_bin_x),end='\r')
      fes[i],der_fes_x[i]=calcFESpoint(grid_cv_x[i])
  else:
    for i in range(grid_bin_x):
      print('   working...  {:.0%} of '.format(i/grid_bin_x),end='\r')
      fes[i]=calcFESpoint(grid_cv_x[i])
else:
  fes=np.zeros((grid_bin_y,grid_bin_x))
  if calc_der:
    der_fes_x=np.zeros((grid_bin_y,grid_bin_x))
    der_fes_y=np.zeros((grid_bin_y,grid_bin_x))
    for i in range(grid_bin_y):
      print('   working...  {:.0%} of '.format(i/grid_bin_y),end='\r')
      for j in range(grid_bin_x):
        fes[i,j],der_fes_x[i,j],der_fes_y[i,j]=calcFESpoint(x[i,j],y[i,j])
  else:
    for i in range(grid_bin_y):
      print('   working...  {:.0%} of '.format(i/grid_bin_y),end='\r')
      for j in range(grid_bin_x):
        fes[i,j]=calcFESpoint(x[i,j],y[i,j])
if args.mintozero:
  fes-=np.amin(fes)
# calculate deltaF
# NB: summing is as accurate as trapz, and logaddexp avoids overflows
if calc_deltaF:
  if not dim2:
    fesA=-kbt*np.logaddexp.reduce(-kbt*fes[grid_cv_x<ts])
    fesB=-kbt*np.logaddexp.reduce(-kbt*fes[grid_cv_x>ts])
  else:
    fesA=-kbt*np.logaddexp.reduce(-kbt*fes[x<ts])
    fesB=-kbt*np.logaddexp.reduce(-kbt*fes[x>ts])
  deltaF=fesB-fesA

### Print to file ###
# backup if necessary
if use_bck:
  cmd=subprocess.Popen('bck.meup.sh -i '+args.outfile,shell=True)
  cmd.wait()
# actual print
f=open(args.outfile,'w')
fields='#! FIELDS '+name_cv_x
if dim2:
  fields+=' '+name_cv_y
fields+=' file.free'
if calc_der:
  fields+=' der_'+name_cv_x
  if dim2:
    fields+=' der_'+name_cv_y
f.write(fields+'\n')
if calc_deltaF:
  f.write('#! SET DeltaF %g\n'%(deltaF))
f.write('#! SET min_'+name_cv_x+' %g\n'%(grid_min_x))
f.write('#! SET max_'+name_cv_x+' %g\n'%(grid_max_x))
f.write('#! SET nbins_'+name_cv_x+' %g\n'%(grid_bin_x))
if period_x==0:
  f.write('#! SET periodic_'+name_cv_x+' false\n')
else:
  f.write('#! SET periodic_'+name_cv_x+' true\n')
if not dim2:
  for i in range(grid_bin_x):
    line=(fmt+'  '+fmt)%(grid_cv_x[i],fes[i])
    if calc_der:
      line+=(' '+fmt)%(der_fes_x[i])
    f.write(line+'\n')
else:
  f.write('#! SET min_'+name_cv_y+' %g\n'%(grid_min_y))
  f.write('#! SET max_'+name_cv_y+' %g\n'%(grid_max_y))
  f.write('#! SET nbins_'+name_cv_y+' %g\n'%(grid_bin_y))
  if period_y==0:
    f.write('#! SET periodic_'+name_cv_y+' false\n')
  else:
    f.write('#! SET periodic_'+name_cv_y+' true\n')
  for i in range(grid_bin_y):
    for j in range(grid_bin_x):
      line=(fmt+' '+fmt+'  '+fmt)%(x[i,j],y[i,j],fes[i,j])
      if calc_der:
        line+=(' '+fmt+' '+fmt)%(der_fes_x[i,j],der_fes_y[i,j])
      f.write(line+'\n')
    f.write('\n')
f.close()
