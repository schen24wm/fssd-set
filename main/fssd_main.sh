#!/bin/bash
## #SBATCH ***
## cd $SLURM_SUBMIT_DIR

# In case something goes wrong at step P and you want to restart from step #Q:
# "rm -r Step{(Q+1)..P}" and "rm -r StepQ/forcecalc" is usually OK.

# Tunable parameters
first_step=1          # The first step (set to 1 for the initial run, higher numbers if restarting)
final_step=100        # Max number of steps
stepsize=0.5          # Step size (unit Bohr)
mixparam=0.36787944   # mixing parameter alpha
# Minimum # of steps before convergence analysis: This parameter depends on geoopt_convergence_analysis.py
minstep=35            # recommend (2*Na + Nb + Nave) or larger
finalcomp=1           # =1: compute a force and energy at the final position. Otherwise: don't compute

# Simulation parameters
noisemag=0.0125       # Gaussian noise to be put onto DFT (unit Ry/Bohr)
sysname=Si            # System name

function F_compute_force() {
  ## Run Dense-K-Grid DFT (Quantum Espresso)
  ## Substitute this part with your own code which generates result.force in the desired format
    mkdir -p QE_DenseKG/
    # QE needs Run_QE.sh, pseudopotential, and POSFILE
    cp  ${basepwd}/qedir_1proc_tmpl/Run_QE.sh QE_DenseKG/
    cp  ${basepwd}/qedir_1proc_tmpl/*.UPF     QE_DenseKG/
    cp  POSFILE_${sysname}                    QE_DenseKG/
    cd  QE_DenseKG/
    ####### entering QE_DenseKG Dir: ######### 
      sh Run_QE.sh >QE.o 2>QE.e
      forcetool__qe -werr > QE.force
      forcetool__scatter QE.force ${noisemag} ../result.force  # to simulate QMC, add a Gaussian-type noise
    ######## leaving QE_DenseKG Dir ########## 
    cd  ..
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mixparam1=$(echo|awk '{printf "%.8f\n",'"${mixparam}"'+1}')
basepwd=$(pwd)
rm -f relpos_all
for istep in `seq ${first_step} ${final_step}` ; do

  oldsteppwd=${basepwd}/Step$((istep-1))
  cd ${oldsteppwd}
  #==============MOVING TO OLD STEP DIRECTORY==============

  #1. Compute force

  mkdir -p forcecalc/
  cp POSFILE_${sysname}-step$((istep-1)) forcecalc/POSFILE_${sysname}
  cd forcecalc/ ; F_compute_force # compute force in forcecalc/ directory
  cd ${oldsteppwd}
  cp forcecalc/result.force step$((istep-1)).force

  #2. Create and move to new directory

  newsteppwd=${basepwd}/Step${istep}
  mkdir -p ${newsteppwd}
  cd ${newsteppwd}
  #==============MOVING TO NEW STEP DIRECTORY==============

  #3. Compute new d (stepPre.force)

  cp ${oldsteppwd}/step$((istep-1)).force .
  if [ ! -f stepPre$((istep-1)).force ] ; then  # If d_(n-1) does not exist, set it to zero [should only apply to n=1]
    forcetool__corr -mult 0 step$((istep-1)).force > stepPre$((istep-1)).force
  else
    cp ${oldsteppwd}/stepPre$((istep-1)).force .
  fi
  forcetool__mixF ${mixparam1} step$((istep-1)).force 1.00000000 stepPre$((istep-1)).force ${mixparam} > stepPre${istep}.force

  #4. Compute new POSFILE

  forcemagn=$(forcetool__norm stepPre${istep}.force)                           # norm of d
  mp=$(echo | awk '{printf "%.16f\n", ('"${stepsize}"')/('"${forcemagn}"');}') # Learning rate (a.k.a. "multiplier")
  cp ${oldsteppwd}/POSFILE_${sysname}-step$((istep-1)) .
  forcetool__applymp POSFILE_${sysname}-step$((istep-1)) stepPre${istep}.force ${mp} POSFILE_${sysname}-step${istep}

  #5. Identify convergence
  #==============MOVING TO ROOT DIRECTORY==============
  cd ${basepwd}
  sh get_single_relpos.sh  Step${istep}/POSFILE_${sysname}-step${istep} > Step${istep}/POSFILE_relpos1_${sysname}-step${istep}
  awk 'NR>9{printf "%22.16f %22.16f %22.16f ",$1,$2,$3}END{print "";}'    Step${istep}/POSFILE_relpos1_${sysname}-step${istep} >> relpos_all
  if [ ${istep} -ge ${minstep} ] ; then
    python geoopt_convergence_analysis.py >gca.out 2>gca.log
    convlvl=$(awk '{print $7}' gca.out)
    cat >&2 gca.out
    if [[ "${convlvl}" == "VLC" ]] ; then
      echo >&2 "Convergence Reached. STOP"
      python geoopt_average_relpos.py > POSFILE_final
      break
    fi
  fi

done

#Compute force and energy at the final position?
if [ ${finalcomp} -eq 1 ] ; then

  cd ${basepwd}/Step${istep}
  mkdir -p forcecalc/
  cp POSFILE_${sysname}-step${istep} forcecalc/POSFILE_${sysname}
  cd forcecalc/ ; F_compute_force # compute force in forcecalc/ directory
  cd ${basepwd}/Step${istep}
  cp forcecalc/result.force step${istep}.force

fi
