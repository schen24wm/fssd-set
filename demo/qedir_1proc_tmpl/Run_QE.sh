#!/bin/bash

#-------------#
# preparation #
#-------------#

#prefix
pref=Si
posfile=POSFILE_${pref}
#fft grid
fftg=( 36 36 36 )  # min 40
#cutoff in the kinetic energy
ecut=25.0
#number of bands (<1/4 * total PWs)
nbnd=40
#K point in the bz
pt=(0.0 0.0 0.0)
#in-run XC
xc=LDA
#still need to manually specify:
#ATOM_SPECIES, K_POINTS cards

lattcnst=`rPu.sh ${posfile} lc`
acount=`rPu.sh ${posfile} ac`
scount=`rPu.sh ${posfile} sc`
#-------------#
# calculation #
#-------------#

time_begin=`date +'%F %T %z'`

#--------------------------------------------------------------------------

# In &system:
suff=SCF-DKG
cat > INPUT-${suff} << ***
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = './',
    outdir='./tmp/',
    prefix='${pref}',
    verbosity='high'
    tprnfor=.true.
 /
 &system
    ibrav=  0,
    celldm(1)=${lattcnst},
    nat= ${acount},
    ntyp= ${scount},
    nbnd= ${nbnd},
    ecutwfc =${ecut},
    input_dft='${xc}',
    occupations='smearing',
    smearing='marzari-vanderbilt',
    degauss=0.002,
    nr1=${fftg[0]},
    nr2=${fftg[1]},
    nr3=${fftg[2]} 
 /
 &electrons
    diago_full_acc = .true.,
    conv_thr = 1.0e-10
!   electron_maxstep = 999
!   diagonalization = 'cg'
!   mixing_beta = 0.10
 /
ATOMIC_SPECIES
 Si  28 14_Si_LDA_25Ry_SRL.UPF  
K_POINTS automatic
6 6 6 0 0 0
ATOMIC_POSITIONS crystal
***
for i in `seq 1 $acount` 
do
  aname=`rPu.sh ${posfile} an $i`
  apos=`rPu.sh ${posfile} ap $i 0`
  echo " ${aname} ${apos}" >> INPUT-${suff}
done
echo CELL_PARAMETERS >> INPUT-${suff}
rPu.sh ${posfile} lv 0 0 >> INPUT-${suff}

mpirun ${QE_DIR}/pw.x < INPUT-${suff} > ${pref}-${suff}.out

#--------------------------------------------------------------------------

time_end=`date +'%F %T %z'`
cat > JOB_SUMMARY << ***
Quantum Espresso Self-Consistent Dense k-point Grid run (QESCDG)
--------------------------------------------------------------------------
Job Prefix: ${pref}
Job begun at: ${time_begin}
    ended at: ${time_end}
***
success_job=`grep 'JOB DONE' ${pref}-${suff}.out | wc -l`
if [ ${success_job} = 1 ]
then
  echo "JOB HAS FINISHED SUCCESSFULLY. ADDITIONAL INFO:" >> JOB_SUMMARY
  printf "MD5 check sum of input file INPUT-${suff}: " >> JOB_SUMMARY
  md5sum INPUT-${suff} | awk '{print $1}' >> JOB_SUMMARY
  printf "FFT Grid guess is ( ${fftg[0]} ${fftg[1]} ${fftg[2]} ) ; QE returns: " >> JOB_SUMMARY
  grep 'Dense  grid' ${pref}-${suff}.out  | awk -F '[\(,\)]' '{printf "( %d %d %d )\n", $2,$3,$4;}' >> JOB_SUMMARY 
  printf "Number of electrons: " >> JOB_SUMMARY
  grep 'number of electrons' ${pref}-${suff}.out | awk '{print $NF}' >> JOB_SUMMARY
  printf "Number of bands guess is ${nbnd} ; QE returns: " >> JOB_SUMMARY
  grep 'Kohn-Sham Wavefunctions' ${pref}-${suff}.out | awk -F '[\(,\)]' '{printf "%d\n", $3}' >> JOB_SUMMARY
  echo "Plane Wave E_cut: ${ecut}" >> JOB_SUMMARY
  printf "Max number of plane waves: " >> JOB_SUMMARY
  grep 'Kohn-Sham Wavefunctions' ${pref}-${suff}.out | awk -F '[\(,\)]' '{printf "%d\n", $2}' >> JOB_SUMMARY
  echo "Exchange Correlation Functional: ${xc}" >> JOB_SUMMARY
  printf "1-electron trial energy: " >> JOB_SUMMARY
  grep 'one-electron contribution' ${pref}-${suff}.out | tail -1 | awk '{print $(NF-1)}' >> JOB_SUMMARY
else
  echo "JOB HAS NOT FINISHED, OR ENCOUNTERED AN ERROR."  >> JOB_SUMMARY
fi

rm -r tmp/
