# functions to use in relation to WDAC 2 (lee paper)
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

def compare_with_domain_orders(query_protein:dict, reference_protein:dict) -> float:
    """
    Compare shared domain pairs between two domain architectures. 
    Order similairty is measured by dividing the shared domains pairs (Qs) by the total domain pairs (Qt).
    Outputs the order similarity (range [0,1])
    """
    Qs = 0
    Qt = 0
    for domain in query_protein:
        if domain in reference_protein:
            Qt += 1
            if query_protein[domain] == reference_protein[domain]:
                Qs += 1
    if (Qt == 0 or Qs == 0):
        return 0
    return Qs / Qt

def compute_weight_score_of_domain(domains:dict, proteins:dict) -> None:
    """
    The domain weight score is the Inverse Abundance Frequency (IAF) multiplied by the Inverse Versatility (IV) of the domain.
    The IAF has the formula IAF(d) - log2(pt/pd) where pt is the number of total proteins and pd is the number of proteins containing domain d.
    The IV has the formula IV(d) = 1 / fd where fd is the number of domains to the left or to the right of domain d, without counting duplicates.

    The input is a dict of domains formatted as {domain1:label1, domain2:label2} where domain is the domain code and label is a string descriptor of the domain,
      and a dict of proteins formatted as {protein1:[domain1, domain2, ...], protein2:[domain1, domain3, ...]} where protein is the protein code and domain1, domain2, ... are the domains in the protein.
    This function writes a file where each line has the domain, the label and its weight score seperated, tab seperated.
    """
    total_proteins = len(proteins)
    domain_count = {}
    domain_families = {}
    for protein in proteins:
        for domain in proteins[protein]:
            if domain not in domain_count:
                domain_count[domain] = 0
            domain_count[domain] += 1
    for protein in proteins:
        for i in range(0, len(proteins[protein])):
            if proteins[protein][i] not in domain_families:
                domain_families[proteins[protein][i]] = set()
            if i != 0:
                domain_families[proteins[protein][i]].add(proteins[protein][i-1])
            if i != len(proteins[protein]) - 1:
                domain_families[proteins[protein][i]].add(proteins[protein][i+1])
    weight_scores = {}
    print(domain_families)
    for domain in domain_count:
        print(domain)
        fd = len(domain_families[domain])
        iaf = np.log2(total_proteins/domain_count[domain])
        iv = 1/fd
        weight_scores[domain] = iaf * iv
    with open("domain_weight_scores.txt", "w") as f:
        for domain in weight_scores:
            f.write(f"{domain}\t{domains[domain]}\t{weight_scores[domain]}\n")



if __name__ == "__main__":
    """
    take in the reference file of protein domain architectures and calculate the weight scores of each domain by calling compute_weight_score_of_domain with a dict containing the domain and its label.
    the format of the reference file is:
    protein_name\tdomain1_name:domain1_label,domain2_name:domain2_label,...
    """
    if len(sys.argv) != 2:
        print("Usage: python formulas.py <reference_file>")
        sys.exit(1)
    reference_file = sys.argv[1]
    domains = {}
    proteins = {}
    with open(reference_file, "r") as f:
        for line in f:
            protein, domains_str = line.strip().split("\t")
            proteins[protein] = []
            for domain_str in domains_str.split(","):
                domain, label = domain_str.split(":")
                proteins[protein].append(domain)
                domains[domain] = label
    compute_weight_score_of_domain(domains, proteins)
