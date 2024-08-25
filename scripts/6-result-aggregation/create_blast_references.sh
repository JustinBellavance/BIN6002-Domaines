while IFS= read -r line; do
	echo "processing $line"
	awk -v myline="$line" '{if (myline==$2){print $5}}' combined_output.txt > reference_per_DIPPA/${line}.txt
done < all_diplonema_in_results.txt
