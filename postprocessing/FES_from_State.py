#! /usr/bin/env python3

### Get the FES estimate used by OPES, from a dumped state file (STATE_WFILE). 1D or 2D only ###
# slightly similar to plumed sum_hills

import sys
import argparse
import numpy as np
import pandas as pd #much faster reading from file
use_bck=False #requires the bck.meup.sh script
#use_bck=True
if use_bck:
  import subprocess

#parser
parser = argparse.ArgumentParser(description='get the FES estimate used by OPES, from a dumped state file (STATE_WFILE). 1D or 2D only')
parser.add_argument('--state',dest='filename',type=str,default='STATE',help='the state file name, with the compressed kernels')
parser.add_argument('--kt',dest='kbt',type=float,required=True,help='the temperature in energy units')
parser.add_argument('--angle1',dest='angle1',action='store_true',default=False,help='the cv1 is an angle in the range [-pi,pi]')
parser.add_argument('--angle2',dest='angle2',action='store_true',default=False,help='the cv2 is an angle in the range [-pi,pi]')
parser.add_argument('--min',dest='grid_min',type=str,required=False,help='lower bounds for the grid')
parser.add_argument('--max',dest='grid_max',type=str,required=False,help='upper bounds for the grid')
parser.add_argument('--bin',dest='grid_bin',type=str,default="100,100",help='number of bins for the grid')
parser.add_argument('--mintozero',dest='mintozero',action='store_true',default=False,help='shift the minimum to zero')
parser.add_argument('--no_der',dest='no_der',action='store_true',default=False,help='skip derivatives to run faster')
parser.add_argument('--all_stored',dest='all_stored',action='store_true',default=False,help='print all the FES stored instead of only the last one')
parser.add_argument('--outfile',dest='outfile',type=str,default='fes.dat',help='name of the output file')
args = parser.parse_args()
#parsing
filename=args.filename
kbt=args.kbt
mintozero=args.mintozero
calc_der=(not args.no_der)
all_stored=args.all_stored
if all_stored:
  if args.outfile.rfind('/')==-1:
    prefix=''
    outfile=args.outfile
  else:
    prefix=args.outfile[:args.outfile.rfind('/')]
    outfile=args.outfile[args.outfile.rfind('/'):]
  if outfile.rfind('.')==-1:
    suffix=''
  else:
    suffix=outfile[outfile.rfind('.'):]
    outfile=outfile[:outfile.rfind('.')]
  outfile=prefix+outfile+'-t%g'+suffix
else:
  outfile=args.outfile
explore='unset'

#get data and check number of stored states
data=pd.read_table(filename,sep='\s+',header=None)
fields_pos=[]
tot_lines=len(data.iloc[:,1])
for i in range(tot_lines):
  if data.iloc[i,1]=='FIELDS':
    fields_pos.append(i)
if len(fields_pos)==0:
  sys.exit(' no FIELDS found in file "'+filename+'"')
if len(fields_pos)>1:
  print(' a total of %d stored states where found'%len(fields_pos))
  if all_stored:
    print('  -> all will be printed')
  else:
    print('  -> only the last one will be printed. use --all_stored to instead print them all')
    fields_pos=[fields_pos[-1]]
fields_pos.append(tot_lines)

for n in range(len(fields_pos)-1):
  print('   working...   0% of {:.0%}'.format(n/(len(fields_pos)-1)),end='\r')
  l=fields_pos[n]
  dim2=False
  if len(data.iloc[l,:])==6:
    cv1name=data.iloc[l,3]
  elif len(data.iloc[l,:])==8:
    dim2=True
    cv1name=data.iloc[l,3]
    cv2name=data.iloc[l,4]
  else:
    sys.exit(' wrong number of FIELDS in file "'+filename+'": only 1 or 2 dimensional bias are supported')
  action=data.iloc[l+1,3]
  if action=="OPES_METAD_state":
    if explore!='no':
      explore='no'
      print(' building free energy from OPES_METAD')
  elif action=="OPES_METAD_EXPLORE_state":
    if explore!='yes':
      explore='yes'
      print(' building free energy from OPES_METAD_EXPLORE')
  else:
    sys.exit(' This script works only with OPES_METAD_state and OPES_METAD_EXPLORE_state')
  if data.iloc[l+2,2]!='biasfactor':
    sys.exit(' biasfactor not found!')
  sf=1 #scaling factor for explore mode
  if explore=='yes':
    sf=float(data.iloc[l+2,3])
  if data.iloc[l+3,2]!='epsilon':
    sys.exit(' epsilon not found!')
  epsilon=float(data.iloc[l+3,3])
  if data.iloc[l+4,2]!='kernel_cutoff':
    sys.exit(' kernel_cutoff not found!')
  cutoff=float(data.iloc[l+4,3])
  val_at_cutoff=np.exp(-0.5*cutoff**2)
  if data.iloc[l+6,2]!='zed':
    sys.exit(' zed not found!')
  Zed=float(data.iloc[l+6,3])
  if explore=='no':
    if data.iloc[l+7,2]!='sum_weights':
      sys.exit(' sum_weights not found!')
    Zed*=float(data.iloc[l+7,3])
  if explore=='yes':
    if data.iloc[l+9,2]!='counter':
      sys.exit(' counter not found!')
    Zed*=float(data.iloc[l+9,3])
  l+=10 #there are always at least 10 header lines
  while data.iloc[l,0]=='#!':
    l+=1
  if l==fields_pos[-1]:
    sys.exit(' missing data!')
#get kernels
  time=float(data.iloc[l,0])
  center_x=np.array(data.iloc[l:fields_pos[n+1],1],dtype=float)
  if dim2:
    center_y=np.array(data.iloc[l:fields_pos[n+1],2],dtype=float)
    sigma_x=np.array(data.iloc[l:fields_pos[n+1],3],dtype=float)
    sigma_y=np.array(data.iloc[l:fields_pos[n+1],4],dtype=float)
    height=np.array(data.iloc[l:fields_pos[n+1],5],dtype=float)
  else:
    sigma_x=np.array(data.iloc[l:fields_pos[n+1],2],dtype=float)
    height=np.array(data.iloc[l:fields_pos[n+1],3],dtype=float)

#set grid
  period_x=0
  grid_bin_x=int(args.grid_bin.split(',')[0])+1
  if args.grid_min is None:
    grid_min_x=min(center_x)
  else:
    grid_min_x=float(args.grid_min.split(',')[0])
  if args.grid_max is None:
    grid_max_x=max(center_x)
  if args.grid_max is None:
    grid_max_x=max(center_x)
  else:
    grid_max_x=float(args.grid_max.split(',')[0])
  if args.angle1:
    if calc_der:
      print(' +++ WARNING: derivatives are not supported for periodic CVs +++')
      calc_der=False
    grid_min_x=-np.pi
    grid_max_x=np.pi
    period_x=2*np.pi
    grid_bin_x-=1
  cv_grid_x=np.linspace(grid_min_x,grid_max_x,grid_bin_x)
  if dim2:
    period_y=0
    if len(args.grid_bin.split(','))!=2:
      sys.exit('two comma separated integers expected after --bin')
    grid_bin_y=int(args.grid_bin.split(',')[1])+1
    if args.grid_min is None:
      grid_min_y=min(center_y)
    else:
      if len(args.grid_min.split(','))!=2:
        sys.exit('two comma separated floats expected after --min')
      grid_min_y=float(args.grid_min.split(',')[1])
    if args.grid_max is None:
      grid_max_y=max(center_y)
    else:
      if len(args.grid_max.split(','))!=2:
        sys.exit('two comma separated floats expected after --max')
      grid_max_y=float(args.grid_max.split(',')[1])
    if args.angle2:
      if calc_der:
        print(' +++ WARNING: derivatives are not supported for periodic CVs +++')
        calc_der=False
      grid_min_y=-np.pi
      grid_max_y=np.pi
      period_y=2*np.pi
      grid_bin_y-=1
    cv_grid_y=np.linspace(grid_min_y,grid_max_y,grid_bin_y)
    x,y=np.meshgrid(cv_grid_x,cv_grid_y)

#calculate
  max_prob=0
  if dim2:
    prob=np.zeros((grid_bin_y,grid_bin_x))
    if calc_der:
      der_prob_x=np.zeros((grid_bin_y,grid_bin_x))
      der_prob_y=np.zeros((grid_bin_y,grid_bin_x))
    for i in range(grid_bin_y):
      print('   working...  {:.0%} of '.format(i/grid_bin_y),end='\r')
      for j in range(grid_bin_x):
        if period_x==0:
          dist_x=(x[i,j]-center_x)/sigma_x
        else:
          dx=np.absolute(x[i,j]-center_x)
          dist_x=np.minimum(dx,period_x-dx)/sigma_x
        if period_y==0:
          dist_y=(y[i,j]-center_y)/sigma_y
        else:
          dy=np.absolute(y[i,j]-center_y)
          dist_y=np.minimum(dy,period_y-dy)/sigma_y
        kernels_ij=height*(np.maximum(np.exp(-0.5*(dist_x**2+dist_y**2))-val_at_cutoff,0))
        prob[i,j]=np.sum(kernels_ij)/Zed+epsilon
        if calc_der:
          der_prob_x[i,j]=np.sum(-dist_x/sigma_x*kernels_ij)/Zed
          der_prob_y[i,j]=np.sum(-dist_y/sigma_y*kernels_ij)/Zed
        if mintozero and prob[i,j]>max_prob:
          max_prob=prob[i,j]
  else:
    prob=np.zeros(grid_bin_x)
    if calc_der:
      der_prob_x=np.zeros(grid_bin_x)
      for i in range(grid_bin_x):
        print('   working...  {:.0%} of '.format(i/grid_bin_x),end='\r')
        if period_x==0:
          dist_x=(cv_grid_x[i]-center_x)/sigma_x
        else:
          dx=np.absolute(cv_grid_x[i]-center_x)
          dist_x=np.minimum(dx,period_x-dx)/sigma_x
        kernels_i=height*(np.maximum(np.exp(-0.5*dist_x*dist_x)-val_at_cutoff,0))
        prob[i]=np.sum(kernels_i)/Zed+epsilon
        if calc_der:
          der_prob_x[i]=np.sum(-dist_x/sigma_x*kernels_i)/Zed
        if mintozero and prob[i]>max_prob:
          max_prob=prob[i]
  if not mintozero:
    max_prob=1

#print out
  if all_stored:
    outfile_n=outfile%time
  else:
    outfile_n=outfile
  if use_bck:
    cmd=subprocess.Popen('bck.meup.sh -i '+outfile_n,shell=True)
    cmd.wait()
  output=open(outfile_n,'w')
  fmt='% 10.6f'
  fields='#! FIELDS '+cv1name
  if dim2:
    fields+=' '+cv2name
  fields+=' file.free'
  if calc_der:
    fields+=' der_'+cv1name
    if dim2:
      fields+=' der_'+cv2name
  print(fields,file=output)
  print('#! SET min_'+cv1name+' '+str(grid_min_x),file=output)
  print('#! SET max_'+cv1name+' '+str(grid_max_x),file=output)
  print('#! SET nbins_'+cv1name+' '+str(grid_bin_x),file=output)
  if period_x==0:
    print('#! SET periodic_'+cv1name+' false',file=output)
  else:
    print('#! SET periodic_'+cv1name+' true',file=output)
  if dim2:
    print('#! SET min_'+cv2name+' '+str(grid_min_y),file=output)
    print('#! SET max_'+cv2name+' '+str(grid_max_y),file=output)
    print('#! SET nbins_'+cv2name+' '+str(grid_bin_y),file=output)
    if period_y==0:
      print('#! SET periodic_'+cv2name+' false',file=output)
    else:
      print('#! SET periodic_'+cv2name+' true',file=output)
    for i in range(grid_bin_y):
      for j in range(grid_bin_x):
        line=(fmt+'  '+fmt+'  '+fmt)%(x[i,j],y[i,j],0-kbt*sf*np.log(prob[i,j]/max_prob))
        if calc_der:
          line+=('  '+fmt+'  '+fmt)%(0-kbt*sf/prob[i,j]*der_prob_x[i,j],0-kbt*sf/prob[i,j]*der_prob_y[i,j])
        print(line,file=output)
      print('',file=output)
  else:
      for i in range(grid_bin_x):
        line=(fmt+'   '+fmt)%(cv_grid_x[i],0-kbt*sf*np.log(prob[i]/max_prob))
        if calc_der:
          line+=('  '+fmt)%(0-kbt*sf/prob[i]*der_prob_x[i])
        print(line,file=output)
  output.close()