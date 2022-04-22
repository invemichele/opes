#! /usr/bin/env python3
#
# script used to generate the data for figS1
#

import numpy as np
import sys

nbins=100
nthetas=100

#Notice that the following is a modified Wolfe-Quapp,
#the original U_wq differs for two coefficients:
#U_wq = x**4 + y**4 - 2 * x**2 - 4 * y**2 + 1 * x * y + 0.3 * x + 0.1 * y
def U(x, y, t=-0.6*np.pi/4):
  if t == 0:
    u = x**4 + y**4 - 2 * x**2 - 4 * y**2 + 2 * x * y + 0.8 * x + 0.1 * y
    return 2 * (u + 9.28)
  else:
    c = np.cos(t)
    s = np.sin(t) 
    return U(c*x-s*y, s*x+c*y, t=0)

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
