import math
import sys

DOMAIN_WEIGHTS = dict()
REF_ARCHITECTURES = dict()
BLASTP_RESULTS = dict()

def dot(X,Y):
	d = 0
	for i in range(len(X)):
		d += (X[i]*Y[i])
	return d


def magnitude(Z):
	m = 0
	for e in Z:
		m += e**2
	return math.sqrt(m)


def similarity(X,Y):
	if magnitude(X) == 0 or magnitude(Y) == 0:
		return 0
	return dot(X,Y) / (magnitude(X)*magnitude(Y))



def order(query_domains:list,reference_domains:list) -> float: # from formulas.py
	Qs = 0
	Qt = 0
	
	sublists1 = [query_domains[i:i+2] for i in range(len(query_domains) - 1)]
	sublists2 = [reference_domains[i:i+2] for i in range(len(reference_domains) - 1)]
	
	for sublist in sublists1:
		if sublist in sublists2:
			Qs += 1
		Qt += 1
	
	for sublist in sublists2:
		if sublist in sublists1:
			Qs += 1
		Qt += 1
	
	return Qs / Qt if Qt != 0 else 0

def load_domains_weights(domain_weights_filepath):
	global DOMAIN_WEIGHTS

	with open(domain_weights_filepath) as f:
		for entry in f:
			entry = entry.split("\t")
			DOMAIN_WEIGHTS[entry[1]] = float(entry[2][:-1]) * 10
			


def load_references(ref_architectures_filepath):
	global REF_ARCHITECTURES

	with open(ref_architectures_filepath) as f:
		for l in f:
			
			architecture_with_id = l[:-1].split("\t")[1].split(",")[:-1]
			architecture = []

			for _domain in architecture_with_id:
				#architecture.append((_domain.split(":")[1]).lower())
				architecture.append(_domain.split(":")[1])
			
			REF_ARCHITECTURES[l[:-1].split("\t")[0]] = architecture
   
def load_blastp_results(blastp_results_filepath):
	global BLASTP_RESULTS

	with open(blastp_results_filepath) as f:
		for l in f:
			l = l.split("\t")
			query = l[0]
			reference = l[1]
			bit_score = float(l[11])
			BLASTP_RESULTS.setdefault(query, []).append({reference: bit_score})


def architecture_to_vector(architecture):
	# return a vector of corresponding domain weights
	# the -1 for unseen/novel domains in our reference db
	return [DOMAIN_WEIGHTS.get(domain, -1) for domain in set(architecture)]


def wdac(input_seqname, input_arch):

	#get unique domains	
	seq_vec = architecture_to_vector(set(input_arch))

	sims = list()
	
	for ref in REF_ARCHITECTURES.keys(): 
		input_arch_set = set(input_arch)
		tmp_vec = seq_vec.copy()

		ref_arch_set = set(REF_ARCHITECTURES[ref])
		ref_vec = architecture_to_vector(REF_ARCHITECTURES[ref])

		combined_set = input_arch_set.union(ref_arch_set)
		combined_list = list(combined_set)
		index_map = {domain: i for i, domain in enumerate(combined_list)}

		new_tmp_vec = [0] * len(combined_list)
		new_ref_vec = [0] * len(combined_list)

		for domain, value in zip(input_arch_set, tmp_vec):
			new_tmp_vec[index_map[domain]] = value

		for domain, value in zip(ref_arch_set, ref_vec):
			new_ref_vec[index_map[domain]] = value

		sim_score = similarity(new_tmp_vec, new_ref_vec)
		order_score = order(input_arch, REF_ARCHITECTURES[ref])
  
		if ((sim_score + order_score) / 2) > 0.75:

			sims.append(
				[
					input_seqname,
					ref,
					sim_score,
					order_score,
					(sim_score + order_score) / 2,
					",".join(input_arch),
					",".join(REF_ARCHITECTURES[ref]),
				]
			)

		#print(ref,s, o, s+o, sep="\t")

	sims.sort(reverse=True, key=lambda entry: entry[4])
	

	# print(sims)
	return(sims)
	# for i in range(50):
	# 	if (sims[i][4] > 1.75):
	# 		print("\t".join(map(str,sims[i])), flush = True)
	



def test():
	input_arch = ["MULE", "OTU"]
	seq_vec = architecture_to_vector(input_arch)

	sims = dict()

	for ref in REF_ARCHITECTURES.keys(): 

		tmp_vec = seq_vec.copy()

		ref_vec = architecture_to_vector(REF_ARCHITECTURES[ref])
	
		ref_vec_len = len(ref_vec)
		tmp_vec_len = len(tmp_vec)

		if ref_vec_len > tmp_vec_len:
			tmp_vec += [0] * (ref_vec_len - tmp_vec_len)
		elif ref_vec_len < tmp_vec_len:
			ref_vec += [0] * (tmp_vec_len - ref_vec_len)
			
		s = similarity(tmp_vec, ref_vec)
		o = order(seq_vec, ref_vec)
		print(ref,s, o, s+o, sep="\t", flush = True)

def getBlastBitScore(query, reference):
	if query in BLASTP_RESULTS:
		for entry in BLASTP_RESULTS[query]:
			if reference in entry:
				return entry[reference]
	return "NA"

if __name__ == '__main__':

	if len(sys.argv) < 4:
		print("Usage: python3 calculate_architecture_similarity.py ref_architectures.tsv domain_weights.tsv input_architectures.tsv > comparison_results.txt")
  		#print("Usage: python3 calculate_architecture_similarity.py ref_architectures.tsv blastp_scores.txt domain_weights.tsv input_architectures.tsv")
		exit()

	print("# Loading domain weights into memory ... ", end="")
	sys.stdout.flush()
	load_domains_weights(sys.argv[2])
	print(" 	done.", flush = True)
 
	# print new way to calculate using blastp
	# print("# Loading BLASTP results into memory ... ", end="", flush=True)
	# load_blastp_results(sys.argv[3])
	# print(" 	done.", flush = True)

	print("# Loading ref architectures into memory ... ", end="")
	sys.stdout.flush()
	load_references(sys.argv[1])
	print(" 	done.", flush = True)

	wdac_results = {}
 
	# i = 0

	with open(sys.argv[3]) as f:
		
		for query in f:
			query = query.split("\t")

			query_seq_name = query[0]
			query_architecture = [d.split(":")[1] for d in query[1][:-1].split(",")]
   			
			sims = wdac(query_seq_name, query_architecture)
			wdac_results[query_seq_name] = sims
	
	for seqname, sims_list in wdac_results.items():
       
		if (len(sims_list) > 0):
			unique_architectures = {}

			#print(sims_list)
			for sims in sims_list:
				unique_architectures.setdefault(sims[6], []).append(sims)
				#print(unique_architectures)
				for unique_architecture, sims_list in unique_architectures.items():
					print("# " + seqname + " 			: " + sims_list[0][5], flush = True)
					print("# Architecture de Référence 	: " + unique_architecture, flush = True)
					print("#prot_id\tmean_score(cos+order/2)", flush=True)
				
					for sims in sims_list:
						print(sims[1], sims[4], flush=True)

						#bitscore = getBlastBitScore(seqname, sims[1])
						#print(sims[1], sims[4], bitscore, flush=True)
					
					# for now, keep it to only one per unique architecture
					break
