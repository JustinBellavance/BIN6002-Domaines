import pandas as pd
import json

#TODO: seperate by species?

# Read the TSV file into a DataFrame
diplonema_architectures = pd.read_csv('../../architectures_Diplonema.tsv', sep='\t', header=None, index_col=False)
uniprot_architectures = pd.read_csv('../../architectures_uniprot.tsv', sep='\t', header=None, index_col=False)


diplonema_architectures.columns = ['key_column', 'value_column']
uniprot_architectures.columns = ['key_column', 'value_column']

print(diplonema_architectures.set_index('key_column').head())
print(uniprot_architectures.set_index('key_column').head())


diplonema_json = diplonema_architectures.set_index('key_column')['value_column'].to_dict()
uniprot_json = uniprot_architectures.set_index('key_column')['value_column'].to_dict()

del diplonema_architectures,uniprot_architectures

def convert_comma_delimited_values_to_list(input_dict):
    updated_dict = {}
    for key, value in input_dict.items():
        if isinstance(value, str) and ',' in value:
            updated_dict[key] = value.split(',')
        else:
            updated_dict[key] = value
    return updated_dict


updated_diplonema = convert_comma_delimited_values_to_list(diplonema_json)
print(diplonema_json)
print(updated_diplonema)

unique_domains= set()
unique_domain_combinations = set()

for key, value in updated_diplonema.items():
    for domains in value:
        unique_domains.add(domains)
    unique_domain_combinations.add(tuple(value))
    
print(unique_domains)
print(unique_domain_combinations)

#updated_refseq = convert_comma_delimited_values_to_list(refseq_json)
#print(refseq_json)
#print(updated_refseq)

updated_uniprot = convert_comma_delimited_values_to_list(uniprot_json)
#print(uniprot_json)
#print(updated_uniprot)

# for every value of set, check if it is present in the values (list) of another dictionary
# if present, increment the counter
# if not present, add it to the dictionary with a counter of 1

# counter_domains_dict_refseq = {}
# counter_domain_combinations_dict_refseq = {}
# for domain in unique_domains:
#     for key, value_list in updated_refseq.items():
#         if domain in value_list:
#             counter_domains_dict_refseq[domain] = counter_domains_dict_refseq.get(domain, 0) + 1
            
# for domain_combination in unique_domain_combinations:
#     for key, value_list in updated_refseq.items():
#         if list(domain_combination) == value_list:
#             counter_domain_combinations_dict_refseq[', '.join(domain_combination)] = counter_domain_combinations_dict_refseq.get(', '.join(domain_combination), 0) + 1
            
# # save results to json
# import json

# with open('counter_domains_dict_refseq.json', 'w') as fp:
#     json.dump(counter_domains_dict_refseq, fp)
    
# with open('counter_domain_combinations_dict_refseq.json', 'w') as fp:
#     json.dump(counter_domain_combinations_dict_refseq, fp)
    
counter_domains_dict_uniprot = {}
counter_domain_combinations_dict_uniprot = {}
for domain in unique_domains:
    for key, value_list in updated_uniprot.items():
        if domain in value_list:
            counter_domains_dict_uniprot[domain] = counter_domains_dict_uniprot.get(domain, 0) + 1
            
for domain_combination in unique_domain_combinations:
    for key, value_list in updated_uniprot.items():
        if list(domain_combination) == value_list:
            counter_domain_combinations_dict_uniprot[', '.join(domain_combination)] = counter_domain_combinations_dict_uniprot.get(', '.join(domain_combination), 0) + 1

with open('counter_domains_dict_uniprot.json', 'w') as fp:
    json.dump(counter_domains_dict_uniprot, fp)
    
with open('counter_domain_combinations_dict_uniprot.json', 'w') as fp:
    json.dump(counter_domain_combinations_dict_uniprot, fp)
