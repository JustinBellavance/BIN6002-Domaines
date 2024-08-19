#!/bin/bash
#SBATCH -J hmmscan_justin
#SBATCH --account=def-gsarah
#SBATCH --time=71:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G

module load hmmer
hmmscan --domtblout hmmscan-results_uniprot.txt --cpu 32 --noali --domE 0.01 pfam_database/Pfam-A.hmm uniprot/uniprot_sprot.fasta