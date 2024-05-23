import numpy as np
import sys
import csv

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

def convert_file_to_dict(file):
    """
    Input: file path to a file containing protein sequences in format: 
    Sequence1_name\tDomain1,Domain2,Domain3
    Sequence2_name\tDomain1,Domain4,Domain5,Domain6

    Convert to dictionary with key = sequence name, value = list of domains.
    """
    dict_of_domains = {}
    with open(file, "r") as f:
        for line in f:
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
    if len(sys.argv) <= 1:
        print("Error: No input file provided.")
        print("Usage: python weight_score_calculation.py refseq_prot_domains.txt")
        sys.exit(1)
    
    file_path = sys.argv[1]
    dict_of_domains = convert_file_to_dict(file_path)
    domains_weight_scores = compute_weight_score_of_domains(dict_of_domains)
   
    # Write the dictionary to a CSV file
    with open('refseq_weight_scores.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for key, value in domains_weight_scores.items():
            csv_writer.writerow([key, value])

    ref_weight_scores = protein_architecture_weight_vector(domain_weight_scores=domains_weight_scores,prot_architectures=dict_of_domains)

    # Write the dictionary to a CSV file
    with open('refseq_prot_weights.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for key, value in ref_weight_scores.items():
            csv_writer.writerow([key, value])

