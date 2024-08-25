#!/bin/bash

# Check if a filename is provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 filename"
    exit 1
fi

# Input file
input_file=$1

# Total number of lines in the file
total_lines=$(wc -l < "$input_file")

# Number of lines per file
lines_per_file=$(( (total_lines + 99) / 100 ))

# Split the file into 500 smaller files without leading zeros in the suffix
split -l "$lines_per_file" --numeric-suffixes=1 --additional-suffix=.txt "$input_file" output_file_

echo "Split the file into 100 files."
