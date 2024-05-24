import numpy as np
import sys

def compute_weight_score_of_domains(dict_of_domains:dict):
    total_proteins = len(dict_of_domains.keys())
    domain_count = {}
    domain_distinct_neighbor = {}

    # Iterating through each domain for each protein sequence
    for prot_domains in dict_of_domains.values():
        prot_domain_count = len(prot_domains)

        # Count number of times domains appear in each sequence
        for index, domain in enumerate(prot_domains):
            if domain not in domain_count:
                domain_count[domain] = 1
            else:
                domain_count[domain] += 1  

            # List for each distinct neighbor of a domain found in sequences
            if domain not in domain_distinct_neighbor:
                domain_distinct_neighbor[domain] = []

            # For domains at the start of the sequence
            if index == 0:
                if prot_domains[index + 1] not in domain_distinct_neighbor[domain]:
                    domain_distinct_neighbor[domain].append(prot_domains[index + 1])
            # For domains between start and end of the sequence
            elif 0 < index < prot_domain_count - 1:
                if prot_domains[index - 1] not in domain_distinct_neighbor[domain]:
                    domain_distinct_neighbor[domain].append(prot_domains[index - 1])
                if prot_domains[index + 1] not in domain_distinct_neighbor[domain]:
                    domain_distinct_neighbor[domain].append(prot_domains[index + 1])
            # For domains at the end of the sequence
            elif index == prot_domain_count - 1:
                if prot_domains[index - 1] not in domain_distinct_neighbor[domain]:
                    domain_distinct_neighbor[domain].append(prot_domains[index - 1])

    # Counting the number of unique neighbors for each domain
    for key, domains in domain_distinct_neighbor.items():
        domain_distinct_neighbor[key] = len(set(domains))

    domain_weight_score = {}
    domain_weight_score.update(domain_count)
    #Calculating IAF
    for key, domains in domain_count.items():
        domain_count[key] = np.log2(total_proteins/domain_count[key]) 
    #Calculating IV 
    for key, domains in domain_distinct_neighbor.items():
        domain_distinct_neighbor[key] = 1/domain_distinct_neighbor[key]

    for key in domain_count.keys() & domain_distinct_neighbor.keys() & domain_weight_score.keys(): 
        domain_weight_score[key] = domain_count[key]*domain_distinct_neighbor[key]

    return(domain_weight_score)

def protein_architecture_weight_vector(domain_weight_scores:dict,prot_architectures:dict):
    prot_architectures_weights = {}
    prot_architectures_weights = {key: [] for key in prot_architectures.keys()}
    
    #Extract weight scores from weight score matrix and append them accordingly for domains present in each protein
    for key, prot_domains in prot_architectures.items(): 
        for domain in prot_domains:
            prot_architectures_weights[key].append(domain_weight_scores[domain])
    return prot_architectures_weights

def compare_domain_architectures(X, Y) -> float:
    """
    Compares the domain architectures of two protein vectors (where each component corresponds to the weight score of a domain)  X and Y.
    Outputs the cosine similarity (range [0,1]) between the two vectors. 1 indicates that the two vectors have the same domains, 0 indicates that they share no domains.
    """
    dot_product = np.dot(X, Y)
    norm_X = np.linalg.norm(X)
    norm_Y = np.linalg.norm(Y)
    return dot_product / (norm_X * norm_Y)



def compare_with_domain_orders(X, Y) -> float:
    """
    Compare shared domain pairs between two domain architectures. 
    In domain evolution, two- or three-domain combinations, called supradomain, are re-used in different protein context, and domain pairs in protein domain architecture occur in only one order, with only about 2% of such pairs occurring on both possible orders.
    Order similairty is measured by dividing the shared domains pairs (Qs) by the total domain pairs (Qt).
    Outputs the order similarity (range [0,1])
    """
    Qs = 0
    Qt = 0
    for i in range(len(X)):
        for j in range(len(Y)):
            if X[i] == Y[j]:
                Qt += 1
                if i == j:
                    Qs += 1
    return Qs / Qt

def convert_file_to_dict(file,remove_spaces=False):
    """
    Input: file path to a file containing protein sequences in format: 
    Sequence1_name\tDomain1,Domain2,Domain3
    Sequence2_name\tDomain1,Domain4,Domain5,Domain6

    Convert to dictionary with key = sequence name, value = list of domains.
    """
    dict_of_domains = {}
    with open(file, "r") as f:

        for line in f:
            if remove_spaces:
                line = line.replace(' ', '')
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split("\t")
            if len(parts) != 2:
                print(f"Skipping malformed line: {line}")
                continue  # Skip lines that don't have exactly two parts
            dict_of_domains[parts[0]] = parts[1].split(",")
    return dict_of_domains


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: No input file provided.")
        print("Usage: python weight_score_calculation.py refseq_architectures.txt query_architectures.txt")
        sys.exit(1)
    
    refseq_arch = convert_file_to_dict(sys.argv[1], remove_spaces=True)
    query_arch = convert_file_to_dict(sys.argv[2], remove_spaces=True)
    dict_of_domains = {**refseq_arch, **query_arch}

    domains_weight_scores = compute_weight_score_of_domains(dict_of_domains)
   
    # Write the dictionary to a .txt file with tab-separated columns
    with open('weight_scores.txt', 'w') as txtfile:
        for key, value in domains_weight_scores.items():
            txtfile.write(f"{key}\t{value}\n")

    refseq_weights = protein_architecture_weight_vector(domain_weight_scores=domains_weight_scores,prot_architectures=refseq_arch)

    # Write the dictionary to a .txt file with tab-separated columns
    with open('refseq_weights.txt', 'w') as txtfile:
        for key, value in refseq_weights.items():
            txtfile.write(f"{key}\t{value}\n")

    query_weights = protein_architecture_weight_vector(domain_weight_scores=domains_weight_scores,prot_architectures=query_arch)

    # Write the dictionary to a .txt file with tab-separated columns
    with open('query_weights.txt', 'w') as txtfile:
        for key, value in query_weights.items():
            txtfile.write(f"{key}\t{value}\n")

    
    #Iterate through query prot sequences
    for query_key, query_weight in query_weights.items():
        for refseq_key, refseq_weight in refseq_weights.items():
            print(query_weight)
            print(refseq_weight)
            print(compare_domain_architectures(query_weight,refseq_weight))
