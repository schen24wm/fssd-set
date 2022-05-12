#!/usr/bin/env python2

import os,sys
import math as M
import numpy as np

def ssqrt(x):
  if x>0:
    return M.sqrt(x)
  elif x>-1.0E-6:   #"eps"
    return 0.0
  else:
    raise Exception('ssqrt takes negative value')

# INPUT FILE can be normal POSCAR or relpos POSCAR

abspos=np.zeros((8,3))


f=open(sys.argv[1],'r')

iline=0
for l in f.readlines():
  if iline == 1:
    lattcnstxy = float(l.split()[0])
  elif iline == 4:
    lattcnstz  = float(l.split()[2])*lattcnstxy
  elif iline >= 8:
    for ii in range(3):
      abspos[iline-8,ii] = float(l.split()[ii])
  iline += 1 
f.close()

# For final processing:
RefMin_relpos=np.zeros((8,3))
RefMin_newpos=np.zeros((8,3))
RefMin_newat=np.zeros(8,dtype=int)
RefMin_dxy    = 10000.0
RefMin_da_z   = 10000.0
RefMin_db_z   = 10000.0
RefMin_d_z    = 10000.0
RefMin_d_all  = 10000.0
RefMin_saCE   = 10000.0
RefMin_saDF   = 10000.0
RefMin_sbCE   = 10000.0
RefMin_sbDF   = 10000.0
RefMin_s_all  = 10000.0

tmp_stdout=sys.stdout ; sys.stdout=sys.stderr #Redirect all print messages from here to stderr
delta_min=10000.0 ; ref_min=-1  # For recognizing RefMin (ref. atom that gives minimum delta)
nerror=0  # Number of errors (if == 8 (number of atoms) => raise Exception)

for refat in range(8):

  relpos=np.zeros((8,3))
  tmppos=np.zeros(3)

  for ii in range(8):
    relpos[ii,:]=(abspos[ii,:]-abspos[refat,:]+0.5)%1.0-0.5  #Result land in [-0.5,+0.5)
  
  newpos=np.zeros((8,3))  # Rearranged position list
  newat=np.zeros(8,dtype=int)
  #H=0,G=1,F=2,E=3,D=4,C=5,B=6,A=7
  
  # Step I: identify the atoms and rearrange
  newpos[0,:]=relpos[refat,:]   #"0" (refat = H = atom 0)
  newat[0]=refat
  
  # Find G (atom 1):
  icount=0
  for ii in range(8):
    if ii != refat:
      if ( relpos[ii,0] > -.25 and relpos[ii,0] < .25 ) and ( relpos[ii,1] > -.25 and relpos[ii,1] < .25 ):
        newpos[1,:]=relpos[ii,:]
        newat[1]=ii
        icount+=1
  if icount != 1:
    print '@@@@WARNING@@@@ Structure warning: cannot identify atom G, for atom H = atom #'+str(refat+1), '; Ref. atom aborted!'
    nerror+=1
    continue
  
  # Find F,E (atom 2,3):
  icount=0
  for ii in range(8):
    if ii != refat:
      if ( relpos[ii,0] < -.25 or  relpos[ii,0] > .25 ) and ( relpos[ii,1] > -.25 and relpos[ii,1] < .25 ):
        newpos[2+icount,:]=relpos[ii,:]
        newat[2+icount]=ii
        icount+=1
  if icount != 2:
    print '@@@@WARNING@@@@ Structure warning: cannot identify atom E,F, for atom H = atom #'+str(refat+1), '; Ref. atom aborted!'
    nerror+=1
    continue
  
  # Find D,C (atom 4,5):
  icount=0
  for ii in range(8):
    if ii != refat:
      if ( relpos[ii,0] > -.25 and relpos[ii,0] < .25 ) and ( relpos[ii,1] < -.25 or  relpos[ii,1] > .25 ):
        newpos[4+icount,:]=relpos[ii,:]
        newat[4+icount]=ii
        icount+=1
  if icount != 2:
    print '@@@@WARNING@@@@ Structure warning: cannot identify atom C,D, for atom H = atom #'+str(refat+1), '; Ref. atom aborted!'
    nerror+=1
    continue
  
  # Find B,A (atom 6,7):
  icount=0
  for ii in range(8):
    if ii != refat:
      if ( relpos[ii,0] < -.25 or  relpos[ii,0] > .25 ) and ( relpos[ii,1] < -.25 or  relpos[ii,1] > .25 ):
        newpos[6+icount,:]=relpos[ii,:]
        newat[6+icount]=ii
        icount+=1
  if icount != 2:
    print '@@@@WARNING@@@@ Structure warning: cannot identify atom A,B (xy directions problem), for atom H = atom #'+str(refat+1), '; Ref. atom aborted!'
    nerror+=1
    continue
  
  # zA <0 and zB >0 must be satisfied. If in the opposite situation, swap A and B
  if newpos[6,2] < 0 and newpos[7,2] > 0:
    tmppos[:]=newpos[7,:]
    newpos[7,:]=newpos[6,:]
    newpos[6,:]=tmppos[:]
    newat[6]+=newat[7]
    newat[7]=newat[6]-newat[7]
    newat[6]-=newat[7]
  if newpos[6,2] <= 0 or newpos[7,2] >= 0:
    print '@@@@WARNING@@@@ Structure warning: cannot identify atom A,B (z direction problem), for atom H = atom #'+str(refat+1), '; Ref. atom aborted!'
    nerror+=1
    continue
  
  #print "                                 Rel. Pos.                                        Abs. Pos."
  for ii in reversed(range(8)):
    #print chr(ord('H')-ii), '(atom #'+str(newat[ii]+1)+')', "  %13.8f  %13.8f  %13.8f" % (newpos[ii,0],newpos[ii,1],newpos[ii,2]), "||", \
                                                            "  %13.8f  %13.8f  %13.8f" % (abspos[newat[ii],0],abspos[newat[ii],1],abspos[newat[ii],2])
  
  # Step II: Calculate xy direction distance d_xy
  
  ideal_xy=np.array([ [0,0], [0,0], [-.5,0], [-.5,0], [0,-.5], [0,-.5], [-.5,-.5], [-.5,-.5] ])
  dx1=0; dx2=0; dy1=0; dy2=0
  for ii in range(8):
    dx1 += ((newpos[ii,0]-ideal_xy[ii,0]+0.5)%1.0-0.5)
    dx2 += ((newpos[ii,0]-ideal_xy[ii,0]+0.5)%1.0-0.5)**2
    dy1 += ((newpos[ii,1]-ideal_xy[ii,1]+0.5)%1.0-0.5)
    dy2 += ((newpos[ii,1]-ideal_xy[ii,1]+0.5)%1.0-0.5)**2
  dxy = ssqrt( dx2 + dy2 - 1.0/8 * dx1**2 - 1.0/8 * dy1**2 )   # see get_relpos_dist.sh. We choose a shift which minimizes the distance.
  #print "d_xy        = ", dxy    # This value should not change for any choice of reference atom
  
  # Step III: Calculate z direction distance (part I = G-A,H-B) da_z  <=== we rearranged d_z^alpha and d_z^beta terms in this alternative calculation
  
  zpos=newpos[:,2]
  
  # Claim: zG-zA and zH-zB is expected to be -1/4 (modulo 1/2)
  da_z = ssqrt(( ( (zpos[0]-zpos[6])%0.5-0.25 )**2 + ( (zpos[1]-zpos[7])%0.5-0.25 )**2 )*0.5)
  # Explaining the formula:   (1) ^ -(-0.25)+0.25=+0.5 here, omitted since %0.5;       (3) ^ norm factor = 1/sqrt(2) from change of coordinates
  #                           (2) -0.25 after %0.5 is to shift the range to [-0.25,0.25)

  #print "d_z^alpha   = ", da_z
  
  # Step IV: Calculate z direction distance (part II = A-B,C-D,E-F,G-H) db_z  <=== we rearranged d_z^alpha and d_z^beta terms in this alternative calculation
 
  # Claim: zA-zB, zC-zD, zE-zF, and zG-zH are expected to be 1/2 (modulo 1) 
  db_z = ssqrt(( ( (zpos[0]-zpos[1])%1.0-0.5 )**2 + ( (zpos[2]-zpos[3])%1.0-0.5 )**2 + ( (zpos[4]-zpos[5])%1.0-0.5 )**2 + ( (zpos[6]-zpos[7])%1.0-0.5 )**2 )*0.5)
  # Explaining the formula:  (1)  ^ -0.5+0.5 here, -0.5 from the expected z-diff                        (2) norm factor = 1/sqrt(2) from change of coordinates ^
  #print "d_z^beta    = ", db_z
  
  d_z = ssqrt (da_z**2 + db_z**2)
  #print "d_z         = ", d_z
  d_all = ssqrt ( (dxy*lattcnstxy)**2 + (d_z*lattcnstz)**2  )
  #print "a_latt(x,y) = ", lattcnstxy
  #print "a_latt(z)   = ", lattcnstz
  #print "d_all[Bohr] = ", d_all
  #print
  
  # Step V: Consider z direction of CE/DF/CF/DE pairs, and the 2D space they construct: this gives sa_<PAIR> and sb_<PAIR>, where <PAIR> = CE/DF/CF/DE
 
  # New (and more consistent) method: use minimum s^2 from CE/DF or CF/DE (Why are we using minimum? because permutations are still permitted in this step)
  sa=[]; sb=[]; s2=[0.0,0.0]
  for ipair in [(3,5),(2,4),(2,5),(3,4)]:
    sa_pair=(zpos[ipair[1]]+zpos[ipair[0]]      +0.25)%0.5-0.25      # zC+zE is expected to be 0   (modulo 1/2)
    sb_pair=(zpos[ipair[1]]-zpos[ipair[0]]-0.25 +0.25)%0.5-0.25      # zC-zE is expected to be 1/4 (modulo 1/2)
    s_pair = ssqrt( sa_pair**2 + sb_pair**2 )
    sa+=[sa_pair,]
    sb+=[sb_pair,]
    if ipair[0]-ipair[1] == -2:
      s2[0]+=s_pair**2  # Pair (C,E) or (D,F)
    else:
      s2[1]+=s_pair**2  # Pair (C,F) or (D,E)
    #print "Pair ("+chr(ord('H')-ipair[1])+","+chr(ord('H')-ipair[0])+") => s^alpha =", "%27.16f"%(sa_pair), ", s^beta =", "%27.16f"%(sb_pair), ", s =", "%27.16f"%(s_pair)

  # In case CF/DE gives a smaller distance, we need to swap E,F; since this is end of processing this will have no effect on anything BUT the summary printing
  if s2[1] < s2[0]:
    newat[2]+=newat[3]
    newat[3]=newat[2]-newat[3]
    newat[2]-=newat[3]
    tmppos[:]=newpos[3,:]
    newpos[3,:]=newpos[2,:]
    newpos[2,:]=tmppos[:]
    saCE=sa[2] ; saDF=sa[3] ; sbCE=sb[2] ; sbDF=sb[3] ; s2min=s2[1]
    #print "Note: E,F swapped"
  else:
    saCE=sa[0] ; saDF=sa[1] ; sbCE=sb[0] ; sbDF=sb[1] ; s2min=s2[0]
 
  s_all  = ssqrt( 0.5 * s2min ) *lattcnstz
  # 1/sqrt(2) from change of coordinates
  #print "s^alpha(CE)  = ", saCE
  #print "s^beta(CE)   = ", sbCE
  #print "s^alpha(DF)  = ", saDF
  #print "s^beta(DF)   = ", sbDF
  #print "s_all[crys.] = ", s_all/lattcnstz
  #print "s_all[Bohr]  = ", s_all
  #print

  # Step VI: Put everything together

  delta_all = ssqrt ( d_all**2 + s_all**2 )
  #print "delta_all    = ", delta_all

  # Check if a new RefMin is found:
  if delta_all < delta_min:
    delta_min=delta_all
    ref_min=refat
    RefMin_relpos[:,:] = relpos[:,:]
    RefMin_newpos[:,:] = newpos[:,:]
    RefMin_newat[:]    = newat [:]
    RefMin_dxy    = dxy
    RefMin_da_z   = da_z
    RefMin_db_z   = db_z
    RefMin_d_z    = d_z
    RefMin_d_all  = d_all
    RefMin_saCE   = saCE
    RefMin_saDF   = saDF
    RefMin_sbCE   = sbCE
    RefMin_sbDF   = sbDF
    RefMin_s_all  = s_all

if nerror == 8:
  print
  print '@@@@ERROR@@@@ Aborted due to structure warnings for every single ref. atom choice.'  
  quit()

sys.stdout=tmp_stdout #Switch back to stdout for output

print   "delta_min    = %27.16f" % delta_min, "<==> Reference atom = #" + str(ref_min+1)
print   "a_latt(xy)   = %27.16f" % lattcnstxy,   ", a_latt(z)    = %27.16f" % lattcnstz
print   "d_all        = %27.16f" % RefMin_d_all
print   "d_xy         = %27.16f" % RefMin_dxy,   ", d_z          = %27.16f" % RefMin_d_z,  ", d_z^alpha    = %27.16f" % RefMin_da_z, ", d_z^beta     = %27.16f" % RefMin_db_z
print   "s_all        = %27.16f" % RefMin_s_all, ", (crys. coord.) %27.16f" % (RefMin_s_all/lattcnstz)
print   "s^alpha(CE)  = %27.16f" % RefMin_saCE,  ", s^beta(CE)   = %27.16f" % RefMin_sbCE, ", s^alpha(DF)  = %27.16f" % RefMin_saDF, ", s^beta(DF)   = %27.16f" % RefMin_sbDF
print "                                 Rel. Pos.                                        Abs. Pos."
for ii in reversed(range(8)):
  print chr(ord('H')-ii), '(atom #'+str(RefMin_newat[ii]+1)+')', "  %13.8f  %13.8f  %13.8f" % (RefMin_newpos[ii,0],RefMin_newpos[ii,1],RefMin_newpos[ii,2]), "||", \
                                                                 "  %13.8f  %13.8f  %13.8f" % (abspos[RefMin_newat[ii],0],abspos[RefMin_newat[ii],1],abspos[RefMin_newat[ii],2])



