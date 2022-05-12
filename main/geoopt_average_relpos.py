#!/usr/bin/env python2
# Note: read relpos_all

import numpy as np
import math as M
import sys

relpos = np.loadtxt('relpos_all')
relpos_ = np.zeros(relpos.shape)
relpos_ave = np.zeros(relpos.shape[1])

try:     # Try finding begin/end step in gca.out
  fgca=open('gca.out','r')
  lgca=fgca.readline().split()
  beg_step=int(lgca[14])+1
  end_step=int(lgca[3])
except IOError:
  if len(sys.argv)>=3:
    pass

if len(sys.argv)>=3:   # Or read these from input parameters (overrides gca.out)
  beg_step=int(sys.argv[1])
  end_step=int(sys.argv[2])  # inclusive. also the reference step
elif len(sys.argv)==2:
  beg_step=int(sys.argv[1])  # Specify begin step only

sys.stderr.write('Average between step %3d ~ %3d\n'%(beg_step,end_step))

for ii in range(beg_step,end_step):

  relpos_[ii,:] = (relpos[ii,:] - relpos[end_step,:]+0.5)%1.0-0.5
  #print ii,relpos_[ii,:]

for ij in range(relpos.shape[1]):

  #print ij,
  relpos_ave[ij] = sum( relpos_[beg_step:end_step+1,ij] )/(end_step+1-beg_step)
  #print relpos_ave[ij],
  relpos_ave[ij] += relpos[end_step,ij]
  #print relpos_ave[ij],
  relpos_ave[ij] = (relpos_ave[ij]+0.5)%1.0-0.5
  #print relpos_ave[ij]

# A comparison with single-steps' results

#refrelpos=np.array([-1.0/3,-1.0/3,-0.081817,-1.0/3,-1.0/3,0.081817])
#for ii in range(beg_step,end_step+1):
#  sumdiff = 0.0
#  for ij in range(relpos.shape[1]):
#    sumdiff += ( (relpos[ii,ij]-refrelpos[ij]+0.5)%1.0-0.5 )**2
#  sys.stderr.write("%3d %22.16f\n" % (ii, M.sqrt(sumdiff)))
#sys.stderr.write("-----------------------\n")
#sumdiff = 0.0
#for ij in range(relpos.shape[1]):
#  sumdiff += ( (relpos_ave[ij]-refrelpos[ij]+0.5)%1.0-0.5 )**2
#sys.stderr.write("AVE %22.16f\n" % M.sqrt(sumdiff))

# Print the average position out like a POSFILE

f=open('Step'+str(end_step)+'/POSFILE_Si-step'+str(end_step),'r')

for ii in range(8):
  sys.stdout.write(f.readline())

ncount=0
for l in f.readlines():
  species = l.split()[3]
  if ncount==0:
    sys.stdout.write("%22.16f %22.16f %22.16f %s\n" % (0,0,0,species))
  else:
    sys.stdout.write("%22.16f %22.16f %22.16f %s\n" % ( tuple(relpos_ave[ncount*3-3:ncount*3]) + (species,) ) )
  ncount+=1 
