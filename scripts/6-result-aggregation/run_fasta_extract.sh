#!/bin/bash
#SBATCH -J split_fasta
#SBATCH --account=def-gsarah
#SBATCH --time=11:58:00
#SBATCH --mem=1G

counter=0
for file in $(ls reference_per_DIPPA); do
    sbatch extract_fasta.sh ../../uniprot_sprot.fasta reference_per_DIPPA/$file split_fastas/$file.faa
    ((counter++))
    if [ $counter -eq 500 ]; then
        wait  # Wait for all jobs to finish before starting the next batch
        counter=0  # Reset the counter
    fi
done

# Wait for any remaining jobs after the loop
wait
