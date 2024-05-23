import sys
import numpy as np
import csv

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
    if len(sys.argv) < 3:
        print("Error: Wrong input provided.")
        print("Usage: python similarity_score_calculation.py refseq_weight_scores.csv refseq_prot_weights.csv query_prot_architectures.txt")
        sys.exit(1)
    
    #Load csv weight scores matrix of DOMAINS
    weight_scores_file = sys.argv[1]
    weight_scores = {}
    # Open the CSV file
    with open(weight_scores_file, mode='r', newline='') as csvfile:
        # Create a reader object with the appropriate delimiter
        csvreader = csv.reader(csvfile, delimiter=';')
        
        # Iterate over each row in the CSV
        for row in csvreader:
            # Use the first column as the key and the second column as the value, converting both to integers
            weight_scores[str(row[0])] = float(row[1])

    #Load query proteins architectures
    query_prot_architectures_file = sys.argv[2]
    query_prot_architectures = convert_file_to_dict(query_prot_architectures_file)
    
    #Give weight scores to query proteins
    query_weight_vector = protein_architecture_weight_vector(query_prot_architectures,weight_scores)

    #Load refseq proteins weight architectures vectors
    refseq_weight_vectors_file = sys.argv[3]
    refseq_weight_vectors = {}
    
    #Load 
    with open(refseq_weight_vectors_file, mode='r', newline='') as csvfile:
        # Create a reader object with the appropriate delimiter
        csvreader = csv.reader(csvfile, delimiter=';')
        
        # Iterate over each row in the CSV
        for row in csvreader:
            # Use the first column as the key and the second column as the value, converting both to integers
            refseq_weight_vectors[str(row[0])] = list(row[1])

    #Iterate through query prot sequences
    for query in query_weight_vector.values():
        for refseq in refseq_weight_vectors.values():
            print(compare_domain_architectures(query,refseq))


    # # Write the dictionary to a CSV file
    # with open('refseq_weight_scores.csv', 'w', newline='') as csvfile:
    #     csv_writer = csv.writer(csvfile)
    #     for key, value in domains_weight_scores.items():
    #         csv_writer.writerow([key, value])