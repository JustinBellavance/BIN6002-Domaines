#!/bin/bash
#SBATCH -J archi_similarity
#SBATCH --account=def-gsarah
#SBATCH --time=02:59:00
#SBATCH --mem=500G

module load python/3.11.5
python calculate_architecture_similarity_uniprot.py ../../architectures_uniprot.tsv ../../domain_weights.tsv ../../architectures_Diplonema.tsv > similarity_results_uniprot_3.0.txt
#python calculate_architecture_similarity.py all_architectures.txt domain_weight_scores.txt blastp/BLAST/blastp_results_merged.txt hp-architectures-hmmscan.tsv > similarity_results_2.0_refseq.txt