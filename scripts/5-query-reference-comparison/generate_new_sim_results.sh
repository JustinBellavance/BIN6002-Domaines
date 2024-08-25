#!/bin/bash
#SBATCH -J archi_similarity
#SBATCH --account=def-gsarah
#SBATCH --time=1:30:00
#SBATCH --mem=5G
#SBATCH --array=1-99

module load python/3.11.5
#python -m pip install numpy
python calculate_architecture_similarity_uniprot.py ../../architectures_uniprot.tsv ../../domain_weights.tsv ../utils/output_file_${SLURM_ARRAY_TASK_ID}.txt > similarity_results_uniprot_${SLURM_ARRAY_TASK_ID}.txt
#python calculate_architecture_similarity.py all_architectures.txt domain_weight_scores.txt blastp/BLAST/blastp_results_merged.txt hp-architectures-hmmscan.tsv > similarity_results_2.0_refseq.txt