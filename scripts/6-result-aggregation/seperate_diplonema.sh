#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_directory>"
    exit 1
fi

# Assign input arguments to variables
INPUT_FASTA=$1
OUTPUT_DIR=$2

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize variables
output_file=""
sequence=""

# Read the FASTA file line by line
while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        # If a sequence header is encountered, save the previous sequence to a file
        if [ -n "$sequence" ]; then
            echo "$sequence" > "$output_file"
        fi
        # Extract the first element of the header (up to the first space or pipe '|')
        header=$(echo "$line" | cut -d' ' -f1 | cut -d'|' -f1 | tr -d '>')
        output_file="$OUTPUT_DIR/${header}.fasta"
        sequence="$line"$'\n'
    else
        # Accumulate sequence lines
        sequence+="$line"$'\n'
    fi
done < "$INPUT_FASTA"

# Save the last sequence to a file
if [ -n "$sequence" ]; then
    echo "$sequence" > "$output_file"
fi

echo "FASTA sequences have been separated into individual files in $OUTPUT_DIR."
