import requests
import os
import sys



# returns the top50 uniprot names for the input Dippa name
def get_uniprot_names_list(dippa_seq_name, similarity_results_filepath): 
	uniprot_list = []
	with open(similarity_results_filepath) as f:
		for line in f:
			if line[0] != "#":
				line_split = line.split("\t")
				if dippa_seq_name in line_split[0]:
					uniprot_list.append(line_split[1])
	return uniprot_list



# Downloads the pdb and saves  
def get_uniprot_pdb(uniprot_name, output_folder): 

	url_base = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb"

	print("Downloading pdb for ",uniprot_name )
	
	r = requests.get(url_base.format(uniprot_name))

	with open(output_folder + "{}.pdb".format(uniprot_name), "wb+") as f:
		f.write(r.content)



if __name__ == "__main__":

	if len(sys.argv) < 3:
		print("usage: python3 get_top50_pdbs.py DIPPA_XXXX similarty_results.csv")
		sys.exit()

	dippa_seq_name = sys.argv[1]
	similarity_results_filepath = sys.argv[2]

	output_folder = dippa_seq_name+"_refs_pdbs/"
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	uniprot_names_list = get_uniprot_names_list(dippa_seq_name, similarity_results_filepath)


	for uniprot_name in uniprot_names_list:
		get_uniprot_pdb(uniprot_name, output_folder)

