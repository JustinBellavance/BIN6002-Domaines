import sys



# This funciton will check if a domain is overlapping with a list
# of domains, using the start and end positions.
def get_overlaps(domain, domains):
	hits = []
	for d in domains:
		if d != domain and ((d[3] >= domain[3] and d[3] <= domain[4]) or (d[4] >= domain[3] and d[4] <= domain[4])):
			hits.append(d)
	return hits



# This function will use the position of the domains to make an architecture
# with overlap handling
# TODO: More sophisticated approach on handling the overlapping domains 
def process_architecture(seq_name ,domains):

	domains = list(sorted(domains, key=lambda x: x[3]))

	if not domains:
		return

	removed_domains = []

	for domain in domains:
		if domain in removed_domains:
			continue

		overlaps = get_overlaps(domain, domains)
		if overlaps:		# overlap handling

			best_evalue_domain = list(sorted(overlaps, key=lambda x: x[5]))[0]

			if best_evalue_domain[5] < domain[5]:
				overlaps.remove(best_evalue_domain)
				removed_domains.append(domain)
				domains.remove(domain)

			for dom in overlaps:
				removed_domains.append(dom)
				domains.remove(dom)

	if len(domains) > 1: 
		#print(seq_name, ",".join(map(lambda domain: domain[1] ,domains)), sep="\t")
		print(seq_name, ",".join(map(lambda domain: domain[1] + ":" + domain[2],domains)), sep="\t")



# This function will parse the cd batch results file and assign the predicted domains
# to each sequence, then submit the domains to process_architecture to get the architecture
def process_cdbatch_results(cdbatch_results_filename): 

	with open(cdbatch_results_filename) as f:

		current_sequence = None	
		current_sequence_domains = []	

		for entry in f:

			if ">" not in entry:
				continue

			entry = entry.split("\t") 

			seq_name = entry[0].split()[2]
			if "pfam" in entry[7]:
				entry[7] = "PF" + entry[7][4:] 

			if not current_sequence: # first run
				current_sequence = seq_name # set focus to seq_name (1st sequence)

			elif current_sequence != seq_name:
				
				process_architecture(current_sequence ,current_sequence_domains) # process architectures
				
				current_sequence = seq_name 	# change the focus to the next sequence 
				current_sequence_domains = []	# reset the domains for the next sequence

			if entry[1] == "specific":
				current_sequence_domains.append((seq_name,entry[7],entry[8],int(entry[3]),int(entry[4]),float(entry[5])))
			elif entry[1] == "superfamily":
				current_sequence_domains.append((seq_name,entry[7],entry[8].split()[0],int(entry[3]),int(entry[4]),float(entry[5])))

		else:
			process_architecture(current_sequence ,current_sequence_domains)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("usage: python", sys.argv[0], "cdbatch_results.tsv")
		exit()
	process_cdbatch_results(sys.argv[1]) 