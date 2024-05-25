# This script serve to calculate domain weights as mentioned in lee article
# ex: python calc-weights.py ref_architectures.tsv
# head ref_architectures.tsv
# seq1	PF0001:dom1,PF0002:dom2
# seq2	PF0002:dom3,PF0003:dom3
# output:
# PF0001	dom1	weight
# PF0002	dom2	weight
# PF0003	dom3	weight

import math
import sys

ARCHITECTURES = []


def IAF(pt, pd):
	return math.log2(pt/pd)


def IV(fd):
	if fd:
		return 1/fd
	return 1 # for solo domains (part of only one mono-domain protein)


def calc_domain_weight(domain):
	# To calculate the weight of a domain, this function will loop through
	# all reference architectures and track the architectures containing the domain
	# and the distinct_neighbors then call IAF and IV

	global ARCHITECTURES

	total_proteins = len(ARCHITECTURES)

	domain_containing_proteins = 0

	distinct_neighbors = set()

	for architecture in ARCHITECTURES:
			
			if domain in architecture:

				domain_containing_proteins += 1

				for i,dom in enumerate(architecture):
					if domain == dom:
						if i != 0: # not the first domain
							distinct_neighbors.add(architecture[ i-1 ]) # add leftside neighbor
						if i != len(architecture) -1 : # not the last domain
							distinct_neighbors.add(architecture[ i+1 ]) # add rightside neighbor
	
	return IAF(total_proteins, domain_containing_proteins) * IV(len(distinct_neighbors)) # * 10 


def get_all_uniq_domains():
	
	domains = set()

	for architecture in ARCHITECTURES:
		for domain in architecture:
			domains.add(domain)

	return domains


def load_architectures(ref_architectures_filepath ,select_name = 1):
	# This function will load all reference architectures into global ARCHITECTURES var
	# We use domain shortnames for architectures presentation but we might use pfam ids
	# by setting select_name = 0

	global ARCHITECTURES

	with open(ref_architectures_filepath) as f:
		for l in f:
			architecture_with_id = l[:-1].split("\t")[1].split(",")[:-1] 	# PF0001:dom1,PF0002:dom2
			# architecture = [_domain.split(":")[select_name] for _domain in architecture_with_id] 	# dom1,dom2
			# architecture = map(lambda dom: dom.split(":")[1], architecture_with_id)

			ARCHITECTURES.append(architecture_with_id)


def main(ref_architectures_filepath):
	load_architectures(ref_architectures_filepath)

	for domain in get_all_uniq_domains():

		print("\t".join(domain.split(":")), calc_domain_weight(domain), sep="\t")


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Input file missing")
		print("usage: python calc-weights.py ref_architectures.tsv")
	else:
		main(sys.argv[1])
	

