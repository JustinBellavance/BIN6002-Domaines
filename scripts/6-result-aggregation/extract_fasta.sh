#!/bin/bash
#SBATCH -J split_fasta
#SBATCH --account=def-gsarah
#SBATCH --time=2:58:00
#SBATCH --mem=2G

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta> <id_list.txt> <output_fasta>"
    exit 1
fi

# Assign input arguments to variables
INPUT_FASTA=$1
ID_LIST=$2
OUTPUT_FASTA=$3

# Create a temporary file to store grep patterns
PATTERN_FILE=$(mktemp)

# Create grep patterns for each ID in the list
while IFS= read -r ID; do
    echo -e ">${ID}\\b"
done < "$ID_LIST" > "$PATTERN_FILE"

# Extract sequences using awk to include all lines until the next header
awk 'BEGIN {RS=">"; FS="\n"} NR==1 {next} {header=$1; seq=""; for (i=2; i<=NF; i++) seq=seq $i} 
    header ~ /('$(
    paste -sd'|' "$PATTERN_FILE"
    )')/ {print ">"header"\n"seq}' "$INPUT_FASTA" > "$OUTPUT_FASTA"

# Remove the temporary pattern file
rm "$PATTERN_FILE"

echo "Subset completed. Extracted sequences are saved in $OUTPUT_FASTA."

