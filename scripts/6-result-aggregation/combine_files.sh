#!/bin/bash

# Output file name
output_file="combined_output.txt"

# Initialize a flag to track whether the header has been added
header_added=false

# Loop through each file matching the pattern
for file in similarity_results_*.txt; do
    if [ "$header_added" = false ]; then
        # Add the entire content of the first file (including the header)
        cat "$file" > "$output_file"
        header_added=true
    else
        # Skip the first line (header) and append the rest to the output file
        tail -n +2 "$file" >> "$output_file"
    fi
done

echo "Files combined into $output_file"

