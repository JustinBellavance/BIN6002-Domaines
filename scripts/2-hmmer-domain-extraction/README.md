
# 2 Using hmmscan

Download the pfam hmm profiles
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
sbatch run_hmmscan_reference.sh
```
