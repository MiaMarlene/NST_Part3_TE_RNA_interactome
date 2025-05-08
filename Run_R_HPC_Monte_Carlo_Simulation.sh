#!/bin/bash
#SBATCH --cpus-per-task 1	# Ensure that all cores are on one ma
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent
#SBATCH --ntasks=10
#SBATCH --array=1-10
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
source /mnt/home3/miska/mb2453/miniforge3/bin/activate part3_env

file="./blood.txt" #The list containing the run accessions
acc=$(awk "NR==$SLURM_ARRAY_TASK_ID" $file)

Rscript MCsim.R ${acc}
