# WDAC pipeline

> [!NOTE]  
> Before executing a command check if the results are already provided to avoid lot of waiting.

#### 1. Prepare db & Calculate domain weights

If not already done gunzip the reference database (download from slack)
```bash
gunzip uniprot_db/uniprot_archs.tsv.gz
```

Than calculate domain weights and save them into a file
```bash
python3 calc-weights.py uniprot_db/uniprot_archs.tsv > uniprot_db/domain_weights.tsv
```

#### 2. Extract hypothetical proteins

```bash
python3 hp_extract.py data/dp-proteome.faa > data/hp-sub.faa
```
Or to include those who dont have domains, add `--all-hp` after the input file. (! might make it default)

```bash
python3 hp_extract.py data/dp-proteome.faa --all-hp > data/hp-all.faa
```

#### 3. Generate domain architectures for input queries
This can be done either by using cd-batch or hmmscan (make sure it is installed)

##### Option 1 : Using cd-batch

Make sure to install request package using
```bash
pip3 install requests
```

Than to submit the sequences to cd-batch run
```bash
python3 cd_batch.py data/hp-sub.faa > data/cd-batch-results.tsv
```

This will output results to `cd-batch-results.tsv` as specified.
To generate architectures from the cd-batch results run
```bash
python3 generate_architectures_from_cdbatch.py data/cd-batch-results.tsv > data/hp-architectures.tsv
```

##### Option 2 : Using hmmscan

Make sure to install hmmer and download the pfam hmm profiles
```bash
sudo apt install hmmer
```
```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz 
mv Pfam-A.hmm.gz data/
```

Run hmmpress on `Pfam-A.hmm`
```bash
hmmpress data/Pfam-A.hmm
```

Then run  hmmrscan on the hypthetical proteins file
```bash
hmmscan --domtblout data/hmmscan-results.tbl --noali -E 0.01 data/Pfam-A.hmm data/hp-sub.faa
```

To generate architectures run the command below
```bash
python3 generate_architectures_from_hmmscan.py data/hmmscan-results.tbl > data/hp-architectures.tsv
```

#### 4. Finaly, Compare architectures

```bash
python3 wdac.py uniprot_db/uniprot_archs.tsv uniprot_db/domain_weights.tsv data/hp-architectures.tsv
```

## Annexe
