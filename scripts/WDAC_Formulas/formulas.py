# functions to use in relation to WDAC (lee paper)
import numpy as np
import sys

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
    if (Qt == 0 or Qs == 0):
        return 0
    return Qs / Qt

def compute_weight_score_of_domain(domains:list, unique_domain_index:dict) -> list:
    """
    The domain weight score is the Inverse Abundance Frequency (IAF) multiplied by the Inverse Versatility (IV) of the domain.
    The IAF has the formula IAF(d) - log2(pt/pd) where pt is the number of total proteins and pd is the number of proteins containing domain d.
    The IV has the formula IV(d) = 1 / fd where fd is the number of distinct domain families adjacent to domain d.
    This function returns a list of the domain weight scores for the given list of domains, within the assigned position of unique_domain_index. Domains that are not present get 0 at that position.
    """
    total_proteins = len(domains)
    domain_count = {}
    domain_families = {}
    for domain in domains:
        if domain not in domain_count:
            domain_count[domain] = 0
        domain_count[domain] += 1
    for domain in domains:
        if domain not in domain_families:
            domain_families[domain] = set()
    for i in range(1, len(domains)):
        domain_families[domains[i]].add(domains[i-1])
        domain_families[domains[i-1]].add(domains[i])
    weight_scores = [0] * len(unique_domain_index)
    print(weight_scores)
    for domain in domain_count:
        fd = len(domain_families[domain])
        iaf = np.log2(total_proteins/domain_count[domain])
        iv = 1/fd
        weight_scores[unique_domain_index[domain]] = iaf * iv
    return weight_scores


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
            line = line.strip().split("\t")
            dict_of_domains[line[0]] = line[1].split(",")
    return dict_of_domains

def create_unique_domain_vector(dict_of_domains:dict) -> dict:
    """
    Create a unique index for each unique domain family inside the dictionary.
    Dictionary format: {domaine1: 0, domaine2: 1, domaine3: 2, ...}
    """
    unique_domain_index = {}
    index = 0
    for seq in dict_of_domains:
        for domain in dict_of_domains[seq]:
            if domain not in unique_domain_index:
                unique_domain_index[domain] = index
                index += 1
    return unique_domain_index

def compare_proteins_to_reference(query:dict, reference:dict):
    """
    Compare a query set of proteins to a reference set of proteins.
    Output: Cosine similarities and order similarities between each query proteins and each reference proteins.
    """
    unique_domain_index = create_unique_domain_vector(reference)
    for seq1 in query:
        for seq2 in reference:
            weight_scores1 = compute_weight_score_of_domain(query[seq1], unique_domain_index)
            weight_scores2 = compute_weight_score_of_domain(reference[seq2], unique_domain_index)
            cosine_similarity = compare_domain_architectures(weight_scores1, weight_scores2)
            order_similarity = compare_with_domain_orders(query[seq1], reference[seq2])
            print(seq1, seq2, cosine_similarity, order_similarity, sep="\t")


def run_all_functions(dict_of_domains:dict):
    """
    Run function on all unique protein combinations in protein_list.
    Don't compare a protein with a protein it has already been compared to. For example, if seq1="A" and seq2="B" have been compared, don't compare seq1="B" with seq2="A".
    Output: File with results of the functions written to it.
    """
    unique_domain_index = create_unique_domain_vector(dict_of_domains)
    items = list(dict_of_domains.items())
    n = len(items)

    for i in range(n-1):
        for j in range(i+1, n):
            (seq1, domains1), (seq2, domains2) = items[i], items[j]
            weight_scores1 = compute_weight_score_of_domain(dict_of_domains[seq1], unique_domain_index)
            weight_scores2 = compute_weight_score_of_domain(dict_of_domains[seq2], unique_domain_index)
            cosine_similarity = compare_domain_architectures(weight_scores1, weight_scores2)
            order_similarity = compare_with_domain_orders(dict_of_domains[seq1], dict_of_domains[seq2])
            print(seq1, seq2, cosine_similarity, order_similarity, sep="\t")



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("input file missing")
        print("ex: python formulas.py domaines.txt > results.txt") # will compare all proteins in the file to each other
        print("OR")
        print("ex: python formulas.py query_domains.txt reference_domains.txt > results.txt") # will compare each protein in query_domains to each protein in reference_domains
    elif len(sys.argv) == 2:
        dict_of_domains = convert_file_to_dict(sys.argv[1])
        print("Sequence1", "Sequence2", "Cosine Similarity", "Order Similarity", sep="\t")
        run_all_functions(dict_of_domains)
    elif (len(sys.argv) == 3):
        query = convert_file_to_dict(sys.argv[1])
        reference = convert_file_to_dict(sys.argv[2])
        print("Sequence1", "Sequence2", "Cosine Similarity", "Order Similarity", sep="\t")
        compare_proteins_to_reference(query, reference)

