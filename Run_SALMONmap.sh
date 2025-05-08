#!/bin/bash
#SBATCH --cpus-per-task 12	# Ensure that all cores are on one ma
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)

#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
source /mnt/home3/miska/mb2453/miniforge3/bin/activate part3_env


#file="/mnt/home3/miska/mb2453/002COMRADES_project/input/COMRADES_lst.txt" #The list containing the run accessions
acc="SRR2147100"
out_dir="/mnt/home3/miska/mb2453/002COMRADES_project/output/SALMON"
#mkdir -p $out_dir/transcripts_quant
index_dir="/mnt/home3/miska/mb2453/001ReferenceMaterial/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/salmonIndex"
input_dir="/mnt/home3/miska/mb2453/002COMRADES_project/input/RNAseq_Senti_2015/cutadapt"

reads1=$input_dir/${acc}_trimmed.fastq.gz
#reads2=$input_dir/${acc}_2_trimmed.fastq.gz

salmon quant -i $index_dir/transcripts_index -l A -p 12 -r $reads1 --validateMappings -o $out_dir/${acc}transcripts_quant

#salmon quant -i $index_dir/transcripts_index -l A -p 12 -1 $reads1 --validateMappings -o $out_dir/${acc}transcripts_quant
