#!/bin/bash
#SBATCH --job-name=r_reg_slid_dFC
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jwk.ng@mail.utoronto.ca
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=12:00:00

#Produce list of subjects to iterate through.
sublist=`cat ${SCRATCH}/FLEXCOG/code/r_sub100.txt`

#Go to submit directory.
cd $SLURM_SUBMIT_DIR

#Get the relevant modules.
module load NiaEnv/2019b
module load intel/2019u4
module load gcc/9.2.0
module load r/3.6.3
module load gnu-parallel/20190322

#Implement the script in parallel on multiple nodes.
parallel 'Rscript --vanilla ${SCRATCH}/FLEXCOG/code/r_slid_dFC.R' {} ::: $sublist ::: 62
