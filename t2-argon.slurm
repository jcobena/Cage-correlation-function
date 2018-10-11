#!/bin/bash
#SBATCH --ntasks=16 --constraint=xeon
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jose.cobena.84@gmail.com


cd "$SLURM_SUBMIT_DIR"

source /usr/usc/intel/16.0/setup.sh
source /usr/usc/openmpi/1.8.8/setup.sh.intel 
mpirun -np 16  /auto/rcf-proj2/jgc/cobenare/lammps/lammps-11Aug17/src/lmp_foo < in-argon-fcc.in

