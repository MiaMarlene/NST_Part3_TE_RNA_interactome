/bin/bash
#SBATCH -n 5                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH --mem=8G                # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mb2453@cam.ac.uk # Email to which notifications will be sent
#SBATCH --array 1-5
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
source /mnt/home3/miska/mb2453/miniforge3/bin/activate

conda activate part3_env

file="./COMRADES_lst.txt" #The list containing the run accessions
acc=$(awk "NR==$SLURM_ARRAY_TASK_ID" $file)
project="COMRADES"

raw="/mnt/home3/miska/mb2453/hpc-work/$project" #The directory containing the fasta files


out="/mnt/home3/miska/mb2453/hpc-work/COMRADES/cutadapt"


cutadapt -q 20 \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-m 30 \
        --output "$out/${acc}_1_trimmed.fastq.gz" \
        --paired-output "$out/${acc}_2_trimmed.fastq.gz" \
        "$raw/${acc}_1.fastq.gz" "$raw/${acc}_2.fastq.gz"







