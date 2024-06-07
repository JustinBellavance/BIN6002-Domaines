#!/bin/bash
#SBATCH -J hmmscan_justin
#SBATCH --account=def-gsarah
#SBATCH --time=24:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --array=2-600%300

module load hmmer
hmmscan --domtblout pfam_database/refseq_domaines_${SLURM_ARRAY_TASK_ID}.txt --noali --cpu 16 -E 0.01 pfam_database/Pfam-A.hmm refseq_protozoa/part-${SLURM_ARRAY_TASK_ID}.faa
