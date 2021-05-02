#!/bin/bash
#SBATCH --job-name=r_vectorizer
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jwk.ng@mail.utoronto.ca
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=12:00:00

#Produce list of subjects to iterate through.
sublist=`cat ${SCRATCH}/FLEXCOG/code/r_sub100.txt` 

#Go to submit directory.
cd $SLURM_SUBMIT_DIR
 
#Get the relevant modules and open virtual environment for libraries.
module load NiaEnv/2019b
module load intel/2019u4
module load python/3.6.8
module load gnu-parallel/20190322
source ${SCRATCH}/FLEXCOG/envs/FlexEnv/bin/activate

#Implement the script in parallel on multiple nodes.
parallel 'python ${SCRATCH}/FLEXCOG/code/r_vectorizer.py' {} ::: $sublist ::: 62 ::: roi_whole.txt
