#!/usr/bin/env python2

import os,sys
import math as M
import numpy as np

def get_dist_vec(lattvec,arr1,arr2): #Give the (arr1-arr2) vector that has minimum |arr1-arr2| when periodic boundary condition is considered.
  l=arr1.shape[0]
  if arr2.shape[0] != l:
    raise Exception('Length mismatch')
  if l%3 != 0:
    raise Exception('Length is not multiple of 3')
  l1=l/3

  distvec=np.zeros(l)

  for ii in range(l1):

    # Get the difference vector for this 3 coordinates
    x0=np.zeros(3)
    for ij in range(3):
      x0[ij] = (arr1[ii*3+ij] - arr2[ii*3+ij]) % 1.0

    # Find the min difference vector length considering periodic condition
    s2min=1.0E+8
    x1vmin=np.zeros(3)
    for ik1 in range(-1,2):
      for ik2 in range(-1,2):
        for ik3 in range(-1,2):
          x1 = x0 + np.array([ik1,ik2,ik3])
          x1v=np.dot(x1,lattvec)
          if sum(x1v**2)<s2min:
            s2min = sum(x1v**2)
            x1vmin = x1v
    distvec[ii*3:ii*3+3]=x1vmin[:]

  return distvec

def conv_to_crys(lattvec,Cart_coords):
  l=Cart_coords.shape[0]
  if l%3 != 0:
    raise Exception('Length is not multiple of 3')
  l1=l/3

  omega = np.linalg.det(lattvec)
  bvec0 = np.cross(lattvec[1,:],lattvec[2,:])/omega
  bvec1 = np.cross(lattvec[2,:],lattvec[0,:])/omega
  bvec2 = np.cross(lattvec[0,:],lattvec[1,:])/omega
  bvec  = np.array([bvec0,bvec1,bvec2])

  Crys_coords=np.zeros(l)
  for ii in range(l1):
    Crys_coords[ii*3:ii*3+3] = np.dot(Cart_coords[ii*3:ii*3+3],np.transpose(bvec))

  return Crys_coords 

X = np.loadtxt('relpos_all')

#======================
lattvec = 8.9875 * np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.104]])
step_uncorrl_period=10  # "Nave". Used for averaging. Larger numbers remove periodic fluctuation.
convergence_criteria=15 # "Nb". Estimate of the minimum "convergence phase" length (so we can confirm a convergence). Larger numbers deal with extremely difficult ravines / skippable local minima.
phase_min_length=5      # "Na". Minimum length for phases (such that we can do accurate statistics for each phase). Must be smaller than q.
# I recommend Na = 4~8, Nave = 2*Na and Nb = 3*Na.  Note that Nave and Nb determines the earliest possible step of detecting a convergence at step M, that is, (M+Nave+Nb).
#======================

try:
  last_step=int(sys.argv[1]) #last step to consider
except IndexError:
  last_step=X.shape[0]-1

vecfin=np.zeros(X.shape[1])
for ii in range(last_step-step_uncorrl_period+1,last_step):
  vecfin += get_dist_vec(lattvec[:,:],X[ii,:],X[last_step,:]) / (step_uncorrl_period-1)
vecfin = conv_to_crys(lattvec,vecfin)
Xfin = X[last_step] + vecfin

new_last_step=last_step-step_uncorrl_period

distlist=np.zeros(new_last_step+1)
fdl=open('relpos_dist_list_Xfin-'+str(step_uncorrl_period),'w')
for ii in range(0,new_last_step+1):
  distlist[ii] = M.sqrt(sum(get_dist_vec(lattvec,X[ii],Xfin)**2))
  sys.stderr.write('%3d %22.16f\n' %(ii,distlist[ii]))
  fdl.write('%3d %22.16f\n' %(ii,distlist[ii]))
fdl.close()
sys.stderr.write('\n')

mincutstep=phase_min_length # recommend set this to p/2
maxcutstep=new_last_step-mincutstep+1

sys.stderr.write('Number of steps in total = %d\n' %(last_step))
sys.stderr.write('Omitted steps at the end (for averaging) p= %d, new step range = 0 ~ %d\n' %(step_uncorrl_period,new_last_step))
sys.stderr.write('Convergence criteria q= %d steps\n' %(convergence_criteria))
sys.stderr.write('Minimum length for phase separation m= %d\n' %(mincutstep))
sys.stderr.write('Cut Step# (CS#) below is NOT included in before-cut, but IS included in after-cut\n')

max_Eratio=0.0
max_step=-1
sys.stderr.write("CS# |      before cut average and error |       after cut average and error          E_bc/E_ac\n")
for cutstep in range(mincutstep,maxcutstep+1):

  s1bc=sum(distlist[:cutstep])
  s2bc=sum(distlist[:cutstep]**2)
  N_bc=cutstep
  A_bc=s1bc/N_bc
  E_bc=M.sqrt( (s2bc/N_bc-A_bc**2) / (N_bc-1) )

  s1ac=sum(distlist[cutstep:])
  s2ac=sum(distlist[cutstep:]**2)
  N_ac=new_last_step+1-cutstep
  A_ac=s1ac/N_ac
  E_ac=M.sqrt( (s2ac/N_ac-A_ac**2) / (N_ac-1) )

  maxdiv=max( abs(distlist[cutstep:]-A_ac)/E_ac )

  sys.stderr.write( '%3d | %16.6f %16.6f | %16.6f %16.6f [ %16.6f ]\n' % (cutstep,A_bc,E_bc,A_ac,E_ac,E_bc/E_ac))  ##,abs(A_bc-A_ac)/E_ac/maxdiv, abs(A_bc-A_ac)/E_bc, abs(A_bc-A_ac)/E_ac, maxdiv)

  if E_bc/E_ac>max_Eratio:
    max_Eratio = E_bc/E_ac
    max_step   = cutstep

sys.stderr.write("\n")
if max_Eratio > 5:
  sys.stderr.write("Convergence: very likely\n")
  conv_lvl='VLC' #very likely convergence
elif max_Eratio > 3:
  sys.stderr.write("Convergence: likely\n")
  conv_lvl='LC' #likely convergence
else:
  conv_lvl='NC' #no convergence
if max_Eratio > 3:
  sys.stderr.write("Max (before cut err) / (after cut err) happens at cut step # %3d .\n" %max_step)
  if max_step >= new_last_step - convergence_criteria +1:
    sys.stderr.write("WARNING: max ratio happens at last %d steps. You might want to wait a few steps to see how it goes.\n" %convergence_criteria)
    conv_lvl+='S' #short
sys.stderr.write("\n")

print 'Analysis at step %3d : ( %s ) Max E_ratio = %16.6f at step %3d ( %3d steps to last step in calculated distance list )' % (last_step, conv_lvl, max_Eratio, max_step, new_last_step+1-max_step)

# recommended criteria for convergence:
# (1) MAX(E_bc/E_ac) > 3 (or 4 or 5)
# (2) this MAX value happens NOT in the last p steps, or the last p/2 steps in display (which only displays thru laststep-p/2-1)
