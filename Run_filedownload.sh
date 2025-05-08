#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent
source /mnt/home3/miska/mb2453/miniforge3/bin/activate

conda activate part3_env


# $1 - Run accession to be downloaded
# $2 - Project ID to name the folder containing Fastq
sample=$1
project="fasta"
mkdir -p /mnt/home3/miska/mb2453/genomes/${project}
out_dir=/mnt/home3/miska/mb2453/genomes/${project}

#out_dir= "$out/${project_name}"
#mkdir -p out_dir
#Choose the first one for PE and second one for SE
#fastq-dump --gzip --split-files --outdir "${out_dir}" "${sample}"
fastq-dump --gzip --outdir "$out_dir" "${sample}"



