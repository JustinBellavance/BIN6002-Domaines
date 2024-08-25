import sys

# returns a dict [uniprot_id] = [goids] 
def get_go_ids(go_ids_filename):
	with open(go_ids_filename) as f:
		next(f)
			
		uniprot_go_ids = dict()

		for line in f:
			line = line[:-1].split("\t")
			uniprot_go_ids[line[0]] = line[1].replace("; ", ",")
	return uniprot_go_ids

if __name__ == '__main__':

	if len(sys.argv) < 3:
		print("usage: python", sys.argv[0], "go_ids.tsv", "comparison_results.tsv")
		sys.exit()

  
	uniprot_go_ids = get_go_ids(sys.argv[1])
  # iterate through the result file and adds go ids at the end 
	with open(sys.argv[2]) as f:
		for line in f:
			if line[0] in ["#", "*"]:
				continue

			uniprot_id = line.split("\t")[3].split("|")[1]
			
			print(line[:-1], uniprot_go_ids[uniprot_id], sep="\t")
