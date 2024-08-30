# WDAC pipeline for D.papillatum hypothetical proteins

#### 1. Extract D.papillatum hypothetical proteins

```bash
subset-only-hypothetical-proteins.py dp-proteome.faa > dippas-hp-sequences.faa
```

#### 2. Domain identification

install hmmer and download the pfam profiles
```bash
sudo apt install hmmer
```
```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz 
```

Run hmmpress on `Pfam-A.hmm`
```bash
hmmpress Pfam-A.hmm
```

Identify domains for the dippas and Uniprot sequences (change input.faa with corresponding file paths)

```bash
hmmscan --domtblout uniprot-hmmscan-results.tbl --noali -domE 0.01 Pfam-A.hmm uniprot.fasta
```
```bash
hmmscan --domtblout dippa-hmmscan-results.tbl --noali -domE 0.01 Pfam-A.hmm dippas-hp-sequences.faa
```

#### 3. Architecture Generation

To be done for both dippas and uniprot sequences

```bash
python3 generate_architectures_from_hmmscan.py uniprot-hmmscan-results.tbl > uniprot_architectures.tsv
```
```bash
python3 generate_architectures_from_hmmscan.py dippa-hmmscan-results.tbl > dippa_architectures.tsv
```

#### 4. Calculate weight scores

Pass Uniprot architectures as input ! 

```bash
python3 calculate-reference-weight-scores.py uniprot_architectures.tsv > domain_weights.tsv
```

#### 5. Finaly, Compare architectures

```bash
python3 comparison_wdac.py
 uniprot_architectures.tsv domain_weights.tsv dippa_architectures.tsv
```
