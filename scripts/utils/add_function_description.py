import sys


def get_uniprot_function_description(go_ids_filename):
	with open(go_ids_filename) as f:
		next(f)
			
		uniprot_function_description = dict()

		for line in f:
			line = line[:-1].split("\t")
			uniprot_function_description[line[0]] = line[1]
	return uniprot_function_description

if __name__ == '__main__':

	if len(sys.argv) < 3:
		print("usage: python", sys.argv[0], "uniprot_function_description.tsv", "comparison_results.tsv")
		sys.exit()
	
	uniprot_function_description = get_uniprot_function_description(sys.argv[1])


	with open(sys.argv[2]) as f:
		next(f)
		for line in f:
			if line[0] in ["#", "*"]:
				continue

			uniprot_id = line.split("\t")[3].split("|")[1]
			
			print(line[:-1], uniprot_function_description[uniprot_id], sep="\t")
