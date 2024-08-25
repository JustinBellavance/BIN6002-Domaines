# Assign input arguments to variables
QUERY_SEQUENCES=$1
REFERENCE_SEQUENCES=$2
OUTPUT_FILE=$3

module load blast+/2.14.1
makeblastdb -in "$REFERENCE_SEQUENCES" -dbtype prot -out tmp/$QUERY_SEQUENCES

# Step 2: Run BLASTP using the custom database
echo "$QUERY_SEQUENCES"
echo "tmp/$QUERY_SEQUENCES"
echo "$REFERENCE_SEQUENCES"
echo "OUTPUT_FILE"
blastp -query "$QUERY_SEQUENCES" -db "tmp/$QUERY_SEQUENCES" -out "$OUTPUT_FILE" -evalue 1000000 -num_threads 4 -max_hsps 1 -outfmt 6

# Optional: Clean up the database files if you don't need them anymore
# Uncomment the line below if you want to remove the database files
# rm custom_blastp_db.*

echo "BLASTP search completed. Results are saved in $OUTPUT_FILE."
