#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task 7
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent

source /mnt/home3/miska/mb2453/miniforge3/bin/activate part3_env

out='/mnt/home3/miska/mb2453/001ReferenceMaterial/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/salmonIndex'

salmon index -p 7 -t $out/gentrome.fa -i $out/transcripts_index --decoys $out/decoys.txt -k 27 --gencode
