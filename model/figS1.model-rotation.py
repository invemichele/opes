#! /usr/bin/env python3
#
# script used to generate the data for figS1
#

import numpy as np
import sys

nbins=100
nthetas=100

def U(x,y,t):
  c=np.cos(t)
  s=np.sin(t)
  return 2*(c*x-s*y)**4+2*(s*x+c*y)**4-4.0*(c*x-s*y)**2-8.0*(s*x+c*y)**2+4*(c*x-s*y)*(s*x+c*y)+1.6*(c*x-s*y)+0.2*(s*x+c*y)

int_grid=np.linspace(-3,3,num=nbins)
deltaF=np.zeros(nthetas)
maximas=np.zeros(nthetas)
thetas_pi4=np.zeros(nthetas)
for n in range(nthetas):
  theta=-np.pi/1.75*n/(nthetas-1)
  thetas_pi4[n]=-theta/np.pi*4
  print('    working... {:.0%}'.format(n/nthetas),end='\r')
  integral=np.zeros(len(int_grid))
  Za=0
  Zb=0
  for i in range(len(int_grid)):
    val=np.array([U(int_grid[i],y,theta) for y in int_grid])
    integral[i]=-np.log(np.trapz(np.exp(-val),int_grid))
    if i<nbins/2:
      Za+=np.exp(-integral[i])
    else:
      Zb+=np.exp(-integral[i])
  integral-=min(integral)
  filename='fes-pi_4_%.3f.data'%(-theta/np.pi*4)
  np.savetxt(filename,np.c_[int_grid,integral],header='x\' fes')
  deltaF[n]=-np.log(Zb/Za)
  center=np.zeros(int(nbins/2))
  it=int(nbins/4)
  for j in range(len(center)):
    center[j]=integral[it]
    it+=1
  maximas[n]=max(center)

np.savetxt('deltaFs.data',np.c_[thetas_pi4,deltaF,maximas],header='theta  deltaF  barrier')
