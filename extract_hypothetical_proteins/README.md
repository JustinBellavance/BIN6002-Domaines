# Example d'utilisation:

## For hypothetical proteins with a domain:
python3 hp-extract.py data/dp-proteome.faa > hp-sub.faa

## For all hypothetical proteins :
python3 hp-extract.py data/dp-proteome.faa --all > hp-sub-all.faa

python3 split4batch.py hp-sub.faa

# To only get the proteins with *.1 in the first column

awk '{if ($1 != "*.1"){print $0}}' architectures.txt > architectures_filtered.txt

# get protozoa protein files

for i in {1..6}
do 
    wget -r https://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/protozoa.${i}.protein.faa.gz
done