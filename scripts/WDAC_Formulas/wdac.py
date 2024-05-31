import math
import sys


DOMAIN_WEIGHTS = dict()
REF_ARCHITECTURES = dict()


def dot_product(X,Y):
	d = 0
	for i in range(len(X)):
		d += X[i]*Y[i]
	return d


def magnitude(Z):
	m = 0
	for e in Z:
		m += e**2
	return math.sqrt(m)


def similarity(X,Y):
	return dot_product(X,Y) / (magnitude(X)*magnitude(Y))


def order(X,Y): # from formulas.py
	Qs = 0
	Qt = 0
	for i in range(len(X)):
		for j in range(len(Y)):
			if X[i] == Y[j]:
				Qt += 1
				if i == j:
					Qs += 1
	if Qt != 0:
		return Qs / Qt
	return 0


def load_domains_weights(domain_weights_filepath):
	#  DOMAIN_WEIGHTS[domain_name] = float: weight

	global DOMAIN_WEIGHTS

	with open(domain_weights_filepath) as f:
		for entry in f:
			entry = entry.split("\t")
			DOMAIN_WEIGHTS[entry[1]] = float(entry[2][:-1]) * 10
			

def load_references(ref_architectures_filepath):
	# Loads reference architectures to a global dict
	# REF_ARCHITECTURES[uniprot_id] = list: architecture

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
	# the -1 is for unseen/novel domains in our reference db
	return [DOMAIN_WEIGHTS.get(domain, -1) for domain in architecture]


def wdac(input_seqname, input_arch):
	# This funciton will take the input architecture and loop 
	# through all ref architectures and calculate similarity and order scores
	# then prints the most top 10 homologous proteins

	seq_vec = architecture_to_vector(input_arch)

	sims = list()
	
	for ref in REF_ARCHITECTURES.keys(): 

		tmp_vec = seq_vec.copy()
		ref_vec = architecture_to_vector(REF_ARCHITECTURES[ref])
	
		ref_vec_len = len(ref_vec)
		tmp_vec_len = len(tmp_vec)

		# adjust the shortest vetcor length to the longest one by adding 0's

		if ref_vec_len > tmp_vec_len: 
			tmp_vec += [0] * (ref_vec_len - tmp_vec_len)
		elif ref_vec_len < tmp_vec_len:
			ref_vec += [0] * (tmp_vec_len - ref_vec_len)
			
		sim_score = similarity(tmp_vec, ref_vec)
		order_score = order(seq_vec, ref_vec)

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

	sims.sort(reverse=True, key=lambda entry: entry[4]) # sort by the (sim_score + order_score) value

	for i in range(10):
		print("\t".join(map(str,sims[i])))


if __name__ == '__main__':

	if len(sys.argv) < 4:
		print("Usage: python3 wdac.py ref_architectures.tsv domain_weights.tsv input_architectures.tsv")
		exit()

	load_domains_weights(sys.argv[2])
	load_references(sys.argv[1])

	with open(sys.argv[3]) as f: 
		for query in f: # loop over each query sequence

			query = query.split("\t")

			query_seq_name = query[0]
			query_architecture = [d for d in query[1][:-1].split(",")]

			wdac(query_seq_name, query_architecture) # check for homologous protiens
