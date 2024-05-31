import math
import sys

DOMAIN_WEIGHTS = dict()
REF_ARCHITECTURES = dict()


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


def sim(X,Y):
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

		
		sim_score = sim(new_tmp_vec, new_ref_vec)
		order_score = order(input_arch, REF_ARCHITECTURES[ref])

		sims.append(
		    [
		        input_seqname,
		        ref,
		        sim_score,
		        order_score,
		        sim_score + order_score,
		        ",".join(input_arch),
		        ",".join(REF_ARCHITECTURES[ref]),
		    ]
		)

		#print(ref,s, o, s+o, sep="\t")

	sims.sort(reverse=True, key=lambda entry: entry[4])


	for i in range(20):
		print("\t".join(map(str,sims[i])))
	



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
			
		s = sim(tmp_vec, ref_vec)
		o = order(seq_vec, ref_vec)

		print(ref,s, o, s+o, sep="\t")



if __name__ == '__main__':

	if len(sys.argv) < 4:
		print("Usage: python3 wdac.py ref_architectures.tsv domain_weights.tsv input_architectures.tsv")
		exit()

	print("Loading domain weights into memory ... ", end="")
	sys.stdout.flush()
	load_domains_weights(sys.argv[2])
	print("done.")

	print("Loading ref architectures into memory ... ", end="")
	sys.stdout.flush()
	load_references(sys.argv[1])
	print("done.")


	with open(sys.argv[3]) as f:
		for query in f:

			query = query.split("\t")

			query_seq_name = query[0]
			query_architecture = [d.split(":")[1] for d in query[1][:-1].split(",")]

			print("-" * 50)
			print("Searching for homologs for ", query_seq_name, query_architecture)
			
			wdac(query_seq_name, query_architecture)
