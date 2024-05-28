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

	#print(seq_name ,domains)
	
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
		#print(seq_name, ",".join(map(lambda domain: domain[1]+ ":" + domain[2] + ":" + str(domain[3]) +"-"+ str(domain[4]),domains)), sep="\t")
		print(seq_name, ",".join(map(lambda domain: domain[1]+ ":" + domain[2],domains)), sep="\t")


# This function will parse the cd batch results file and assign the predicted domains
# to each sequence, then submit the domains to process_architecture to get the architecture
def process_hmmscan_results(hmmscan_results_filename): 

	hmm_all_entries = []

	with open(hmmscan_results_filename) as f:

		current_sequence = None	
		current_sequence_domains = []	

		for entry in f:

			if "#" == entry[0]:
				continue
			entry = entry.split()
			try:
				hmm_all_entries.append([entry[i] if i != 1 else entry[i].split(".")[0] for i in  [0, 1, 3, 6, 17, 18]])
			except:
				pass

		hmm_all_entries.sort(key = lambda e:e[2])
		
		for entry in hmm_all_entries:

			seq_name = entry[2]
			#print(seq_name)
			if not current_sequence: # first run
				current_sequence = seq_name # set focus to seq_name (1st sequence)

			elif current_sequence != seq_name:
				
				process_architecture(current_sequence ,current_sequence_domains) # process architectures
				
				current_sequence = seq_name 	# change the focus to the next sequence 
				current_sequence_domains = []	# reset the domains for the next sequence
				# seqname pfamid pfam name start end evalu
			current_sequence_domains.append((seq_name,entry[1],entry[0],int(entry[4]),int(entry[5]),float(entry[3])))

		else:
			process_architecture(current_sequence ,current_sequence_domains)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("usage: python", sys.argv[0], "hmmscan.results.tbl")
		exit()
	process_hmmscan_results(sys.argv[1]) 