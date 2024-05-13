# Nous allons commencer avec selectionner que les séquences hypothetique en utilsant le script extract_hypothetical_proteins/hp-extract.py

python3 ../extract_hypothetical_proteins/hp-extract.py ../dp-proteome.faa --all > ../hp-sub-all.faa

# Nous allons utiliser SMART normal pour les domaines (https://smart.embl.de/smart/set_mode.cgi?NORMAL=1#)
## Selon eux (https://smart.embl.de/help/FAQ.shtml), il n'a pas de bonne facon de faire du batch serving pour sequences inconnus, donc il faut utilsé ce script perl. (Perl est probalement déja installé sur systèmes linux). Je recommande bioperl pour les modules utilisés dans ce script (

wget https://smart.embl.de/help/SMART_batch.pl
sudo apt-get install bioperl bioperl-run

./SMART_batch.pl --help
./SMART_batch.pl --includePfam --inputFile ../hp-sub.faa --outputDirectory validation_test #pour tester
./SMART_batch.pl --includePfam --inputFile ../hp-sub-all.faa --outputDirectory validation_test #actuel

# Avec les résultats du script, il faut un script python pour bien selectionner les domaines avec des bons e-values. Ils ont un different status que les autres. (moins que 0.01)
