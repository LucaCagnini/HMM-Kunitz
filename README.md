
# Profile HMM Pipeline for Kunitz-Type Protease Inhibitor Domain (Pfam: PF00014)

This repository implements a **computational pipeline** for building and evaluating a Profile Hidden Markov Model (HMM) designed to detect the **Kunitz-type protease inhibitor domain** (Pfam ID: PF00014) within protein sequences.

The project was developed as part of the _Laboratory of Bioinformatics 1_ course (MSc in Bioinformatics, University of Bologna) and integrates **structural bioinformatics**, **sequence analysis**, and **statistical model evaluation**.

---
## Table of Contents

- [Objectives](#objectives)
- [Requirements](#requirements)
- [Pipeline Execution](#pipeline-execution)
  - [1. Extract Kunitz Sequences](#1-extract-kunitz-domain-sequences-from-uniprot)
  - [2. Perform MSA](#2-perform-msa)
  - [3. Build the HMM](#3-Build-the-hmm)
  - [4. Search with the Model](#4-search-with-the-model)
  - [5. Format Results](#5format-results-into-class-files)
  - [6. Evaluate the Model](#6-evaluate-the-model)
- [Output](#-output)
- [Author](#-author)
- [License](#-license)

  ---

## objectives

- Build a structure-based multiple sequence alignment of known Kunitz-domain proteins.
- Filter redundancy using BLAST and CD-HIT.
- Create and calibrate an HMM model using that alignment.
- Evaluate the model’s performance using curated positive (true Kunitz) and negative (non-Kunitz) datasets.
- Output metrics such as confusion matrix, accuracy, sensitivity, and specificity. 

---

## Requirements

All the requirement packages needs to be installed before. This can be done in a dedicated conda enviroment.

```bash
conda create -n kunitz_env python=3.10
conda activate kunitz_env
conda install -c bioconda
conda install -c hmmer
conda install -c blast
conda install -c biopython
conda install cd-hit
```

## Pipeline Execution

### 1. *Extract Kunitz domain sequences from UniProt:*

From the advance search function in UniProt it is possible to download the sequences of all the proteins showing a Kunitz domain, with a specific query:
   
   ° Pfam ID: PF 00014
   
   ° Resolution ≤ 3.5 Angstrom
   
   ° 45 ≤ length ≥ 80
   
Query:

```
Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
```

Then we can extract the IDs from the results file (rcsb_pdb_custom_report20250410062557.csv) and create a fasta file:

```
cat rcsb_pdb_custom_report20250410062557.csv| tr -d '"' | awk -F ',' '{if (length($2)>0) {name=$2}; print name ,$3,$4,$5}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz_customreported.fasta
```

Then, we cluster the sequences using CD-HIT, at 90% identity threshold.
CD-HIT is a greedy incremental algorithm that starts with the longest input sequence as the first cluster representative and then processes the remaining sequences from long to short to classify each sequence as a redundant or representative sequence based on its similarities to the existing representatives (Fu et al., 2012).


```
cd-hit -i pdb_kunitz_customreported.fasta -o pdb_kunitz_customreported.clstr -c 0.9

```
We then take the file an analyse wether we can need to cut clusters away, before proceeding with our Multiple Sequence Alignment. results file: pdb_kunitz_customreported_trimmed.fasta 
We then extract the most representative ID from each cluster

```

clstr2txt.pl pdb_kunitz_customreported_filtered.clstr > pdb_kunitz.clusters.txt

```
we then take the ids from that file and we store them in a FASTA file

```

awk '$5 == 1 {print $1}' pdb_kunitz.clusters.txt > pdb_kunitz_rp.ids

```

and prepare a file for the MSA

```
for i in $(cat pdb_kunitz_rp.ids); do   grep -A 1 "^>$i" pdb_kunitz_customreported.fasta | head -n 2 >> pdb_kunitz_rp.fasta; done

```

pdb_kunitz_rp.fasta contains only the represantative sequences selected from CD-HIT

### 2. Perform MSA

We used the PDBeFold Multi alignment Tool
[link text] (https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver)

we download the PDEBefold alignment as .ali file and then format it

```
awk '{if (substr($1,1,1)==">") {print "\n" toupper($1)} else {printf "%s", toupper($1)}}'pdb_kunitz_rp.ali >pdb_kunitz_rp_formatted.ali
```

### 3. Bulid the HMM

We activate our dedicated conda enviroment and build the HMM model using the output file from the previous MSA

```
conda activate 
hmmbuild structural_model.hmm pdb_kunitz_rp_formatted.ali

```

### 4. Search with the model 

We run hmmsearch on all positive and negative FASTA files and create a tabular output. We generate .out files that contains the E-values computed on the HMM of the sequence alignment.

```
hmmsearch -Z 1000 --max --tblout pos_1_seqali.out pdb_kunitz_rp_formatted.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2_seqali.out pdb_kunitz_rp_formatted.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1_seqali.out pdb_kunitz_rp_formatted.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2_seqali.out pdb_kunitz_rp_formatted.hmm neg_2.fasta
```

we create a database with the sequences, downloading them from PDB:

-all_kunitz.fasta

-human_kunitz.fasta

-human_not_kunitz.fasta

It is necessary to remove redundant proteins from our tests set. 

```
makeblastdb -in all_kunitz_uniprot.fasta -dbtype prot -out all_kunitz_uniprot.fasta

blastp -query pdb_kunitz_rp.fasta -db all_kunitz_uniprot.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
```
we remove all the proteins with high identity to the one with witch we have created our HMM model. 

```
grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > to_remove.ids

grep ">" all_kunitz_uniprot.fasta | cut -d "|" -f 2 > all_kunitz.id

comm -23 <(sort all_kunitz.id) <(sort to_remove.ids) >to_keep.ids

```

### 5. Format results

we need to obtain the sequences of the proteins we keep, to do so we use the script *get_seq.py*
and we need to divide in negative and positive set. 

```
python3 get_seq.py to_keep.ids all_kunitz.fasta ok_kunitz.fasta

grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 >sp.id

comm -23 <(sort sp.id) <(sort all_kunitz.id) >sp_negs.ids

python3 get_seq.py sp_negs.ids uniprot_sprot.fasta sp_negs.fasta

```
we need to randomise the positive and negative set

```
sort -R sp_negs.ids > random_sp_negs.ids
sort -R to_keep.ids > random_ok_kunitz.ids

head -n 184 random_ok_kunitz.ids > pos_1.ids
tail -n 184 random_ok_kunitz.ids > pos_2.ids

head -n 286417 random_sp_negs.ids >neg_1.ids
tail -n 286417 random_sp_negs.ids >neg_2.ids

```

### 6.Evaluate the model

We use hmmsearch ad use as input our model and our positives/negatives sets. 

```
hmmsearch -Z 1000 --max --tblout pos_1.out structural_model.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out structural_model.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out structural_model.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out structural_model.hmm neg_2.fasta

```
We then convert the output in a classification format 

```

grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_2.class

```
we create two set from the previous files

```
cat pos_2.class neg_2_hits.class >set_2.class
cat pos_2.class neg_2_hits.class > set_2.class
```
we use the file performance.py to value the performance

```

python3 performance.py set_2.class 1e-5 >results_set_1.txt
python3 performance.py set_2.class 1e-5 >results_set_2.txt

```

we can run the command, using different thresholds for the E value, evalutaing the changes of our performances. 

```

for i in $(seq 1 10); do   python3 performance.py set_1.class 1e-$i; done | sort -nrk 6 > diff_threshold_set1.txt
for i in $(seq 1 10); do   python3 performance.py set_2.class 1e-$i; done | sort -nrk 6 > diff_threshold_set2.txt

```

 



