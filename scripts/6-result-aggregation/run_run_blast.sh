#!/bin/bash
#SBATCH -J blastp
#SBATCH --account=def-gsarah
#SBATCH --time=00:58:00
#SBATCH --mem=1G
#SBATCH --ntasks=1                  
##SBATCH --cpus-per-task=4           
#SBATCH --array=1-5440%500

query_files=($(ls -1 diplonema_faa))
reference_files=($(ls -1 split_fastas))

query_file=${query_files[$((SLURM_ARRAY_TASK_ID-1))]}
reference_file=${reference_files[$((SLURM_ARRAY_TASK_ID-1))]}

query_file="${reference_file%.txt.faa}.fasta"

output_file="$query_file"

echo "$query_file" 
echo "$reference_file"
echo "$output_file"



bash run_blast.sh diplonema_faa/$query_file split_fastas/$reference_file BLASTP_output/$output_file
