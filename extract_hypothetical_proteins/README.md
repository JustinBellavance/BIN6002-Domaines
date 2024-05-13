# Example d'utilisation:

## For hypothetical proteins with a domain:
python3 hp-extract.py data/dp-proteome.faa > hp-sub.faa

## For all hypothetical proteins :
python3 hp-extract.py data/dp-proteome.faa --all > hp-sub-all.faa

python3 split4batch.py hp-sub.faa