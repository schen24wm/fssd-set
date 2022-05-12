This repository contains scripts that implements the `FSSDxSET` algorithm, which aims to realize fast and robust structural optimization of solids with accurate but noisy and expensive (e.g. many-body) forces.

A preprint of the paper that describes the algorithm can be found on [arXiv](https://arxiv.org/abs/2204.12074) and is submitted to Nature Computational Science.

These scripts are not standalone in themselves. A force computing software is required to generate forces. An additional script may be required to process the output to the specified format.

# Software Requirements

The scripts are tested on CentOS Linux 7.

Software requirements include:

- GNU Bash > 3.0 (Tested on 4.2.46)
- Python2 > 2.7 (Tested on 2.7.13)
- External software for force computation (e.g. Quantum Espresso, PW-AFQMC)

Python2 scripts use `numpy` as a dependency.

If using Quantum Espresso to simulate many-body forces (as in the demo),
- Quantum Espresso > 5.0.1 (Tested on 5.0.1, 6.3 and 6.8)

# Installation

No actual installation is needed if the software requirement is met, but a few setup of the running environment are needed.

- Scripts in `bin-addon` needs to be added to a directory in `PATH` and be made executable.
- Set the environment variable `SYSTMP` to an empty directory. This directory will store some (hidden) temporary files of the bin-addon/ scripts.

The entire process can be done within 5 minutes.

# Demo

In the `demo` directory, we provide a two-stage example of `FSSDxSET` in silicon, based on the diamond to beta-tin phase transition example in our paper.

## Preparation

Before running the demo, you should copy everything in `main` to a new empty directory (e.g. `project`).

You also need to copy the initial position file (`demo/POSFILE_Si-step0`) to `project/Step0`.

To use Quantum Espresso, copy the directory `demo/qedir_1proc_tmpl` to `project` as well.

Tweak running parameters in `fssd_main.sh` and `geoopt_convergence_analysis.py` if needed.

## Running

To run the demo (first stage), in the `project` directory, run
```
sh fssd_main.sh
```
However, running with MPI and multiple processors is recommended, as the external force computations (which are called by the script) are usually computationally costly.

## If a second stage is desired

Copy the scripts in `main` to a directory for the second stage (e.g. `Stage2`). Copy `POSFILE_final` in stage 1 to e.g. `Stage2/Step0`.

Rerun `fssd_main.sh` at the `Stage2` directory. Note that following the SET approach, you need to reduce the step size and noise size. For the latter, if you are doing a simulation, you can reduce the scale of the add-on Gaussian noise. If you are running with an actual noisy force program (e.g. PW-AFQMC), you need to change the amount of samples or computational time, so that the target error bar (computed from central limit theorem) reduce by the desired amount.

## Expected behavior

In stage 1, with the preset parameters, the convergence should be identified around (`N` in Appendix E of the paper) step 55 and should be reached around (`m` in Appendix E of the paper) step 25\~30.

The `POSFILE_final` in stage 1 should match with `demo/POSFILE_Si-betatin` with very small fluctuations (max crystal coordinate deviation \~ 0.01).

If a stage 2 is run, this scale of fluctuation should be smaller than stage 1 (depending on the parameters chosen).

## Expected running time

The running time (excluding the call time of the force code, e.g. Quantum Espresso) is about 6 minutes for a stage (with 50~60 steps).
The corresponding Quantum Espresso run needs about 3 hours when run serially.

# Instructions for use

## Running with your own solid

To do optimizations with your own solid, make an initial position file with the same format and replace `POSFILE_Si-step0`. Please note that this position file is almost identical to VASP's POSCAR file, except that the lattice constant is in Bohr, not Angstrom. The atom positions are in crystal coordinates.

Check all scripts in `project` and `project/qedir_1proc_tmpl`, to remove dependencies on the solid system (e.g. use of Si pseudopotential and lattice constants).

## Running with your own force code

To run with your own force code, pay attention to the force output format of the demo, e.g. `demo/step0_demo.force`. You need to format the force output of your force code to the same format to make the scripts usable.

If `Nat` atoms are present, the force file contains `Nat` lines, with each line being something like
```
     atom    <N> type  <T>   force =    <Fx>    <Fy>   <Fz>   ,err =    <Ex>    <Ey>    <Ez>
```
where `<N>` is the atom number (beginning with 1), `<T>` is the atom type number (beginning with 1). `<Fx> <Fy> <Fz>` is the force acting on the atom in _Cartesian_ coordinates and _Ry/Bohr_ units. `<Ex> <Ey> <Ez>` contains the corresponding error bars for the forces.
