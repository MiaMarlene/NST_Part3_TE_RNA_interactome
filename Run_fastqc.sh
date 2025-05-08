#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent
source /mnt/home3/miska/mb2453/miniforge3/bin/activate

conda activate part3_env

out=$1


fastqc --out-dir "${1}/fastqc/" ${1}/*.fastq.gz
