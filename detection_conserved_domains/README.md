# Nous allons commencer avec selectionner que les séquences hypothetique en utilsant le script extract_hypothetical_proteins/hp-extract.py

python3 ../extract_hypothetical_proteins/hp-extract.py ../dp-proteome.faa --all > ../hp-sub-all.faa

num_batches=$(grep ">")

# Nous allons utiliser SMART normal pour les domaines (https://smart.embl.de/smart/set_mode.cgi?NORMAL=1#)
## Selon eux (https://smart.embl.de/help/FAQ.shtml), il n'a pas de bonne facon de faire du batch serving pour sequences inconnus, donc il faut utilsé ce script perl. (Perl est probalement déja installé sur systèmes linux). Je recommande bioperl pour les modules utilisés dans ce script

wget https://smart.embl.de/help/SMART_batch.pl
sudo apt-get install bioperl bioperl-run

./SMART_batch.pl --help
./SMART_batch.pl --includePfam --inputFile ../hp-sub.faa --outputDirectory validation_test #pour tester
./SMART_batch.pl --includePfam --inputFile ../hp-sub-all.faa --outputDirectory new_domains #actuel

# Mais c'est trop lent pour des queries plus que 100.

# Avec  CD-batch, un autre perl script, il faut changer le format du fichier un peu. (une sequence par ligne, pas de header)

./bwrpsb.pl < ../hp-sub.faa

# Avec les résultats du script, il faut un script python pour bien selectionner les domaines avec des bons e-values. Ils ont un different status que les autres. (moins que 0.01)

rm multi_domain_query.txt
for i in {0..315}
do
    grep -v "^#" part-${i}_hitdata.txt| awk 'NR>2{print $1}' | sort | uniq -c | awk '{if ($1 > 1){print $0}}' >> multi_domain_query.txt
done

# pour combiner les resultats de domaines 
rm domain_results_combined.txt
head -n 8 part-0_hitdata.txt > domain_results_combined.txt
for i in {0..1178}
do
    grep -v "^#" part-${i}_hitdata.txt| awk 'NR>2{print $0}' >> domain_results_combined.txt
done


# Pour telecharger les FASTA de PFAM:

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz

# run all files

for i in {0..333}
do
    python cd_batch.py ../refseq_protozoa/part-${i}.faa > ../Batch-CD-protozoa_results/part-${i}_hitdata.txt
done

for i in {334..667}
do
    python cd_batch.py ../refseq_protozoa/part-${i}.faa > ../Batch-CD-protozoa_results/part-${i}_hitdata.txt
done

for i in {667..1000}
do
    python cd_batch.py ../refseq_protozoa/part-${i}.faa > ../Batch-CD-protozoa_results/part-${i}_hitdata.txt
done

for i in {1000..1178}
do
    python cd_batch.py ../refseq_protozoa/part-${i}.faa > ../Batch-CD-protozoa_results/part-${i}_hitdata.txt
done