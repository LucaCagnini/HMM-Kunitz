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
  - [3. Build the HMM](#3-build-the-hmm)
  - [4. Model Testing](#4-model-testing)
  - [5. Format Results](#5-format-results)
  - [6. Model Evaluation](#6-model-evaluation)
- [Output](#output)
- [Author](#author)
- [License](#license)

---

## Objectives

- Build a structure-based multiple sequence alignment of known Kunitz-domain proteins.
- Filter redundancy using BLAST and CD-HIT.
- Create and calibrate an HMM model using that alignment.
- Evaluate the model’s performance using curated positive (true Kunitz) and negative (non-Kunitz) datasets.
- Apply output metrics such as confusion matrix, accuracy, sensitivity, and specificity. 

---

## Requirements

Before running the pipeline, all required programmes were installed in a Conda environment:

```bash
conda create -n kunitz_env python=3.10
conda activate kunitz_env
conda install -c bioconda hmmer blast biopython
conda install -c conda-forge cd-hit
```

to create and visualise plot, the following python libraries have been used. 

```bash
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import nunpy as np
from graphviz import Digraph
```


## Pipeline Execution

N.B. in this repositry have been uploaded all the output files of the following script, as well as the python programmes used under the  files/ folder. 

### 1. *Extract Kunitz domain sequences from UniProt:*

From the advance search function in UniProt it is possible to download the sequences of all the proteins showing a Kunitz domain, with a specific query:
   
   ° Pfam ID: PF 00014
   
   ° Resolution ≤ 3.5 Angstrom
   
   ° 45 ≤ length ≥ 80
   
Query:

```
Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
```

Then the IDs from the results file (rcsb_pdb_custom_report20250410062557.csv) are extracted and a fasta file is created:

```
cat rcsb_pdb_custom_report20250410062557.csv| tr -d '"' | awk -F ',' '{if (length($2)>0) {name=$2}; print name ,$3,$4,$5}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz_customreported.fasta
```

The sequences are clustered using CD-HIT, at 90% identity threshold.
CD-HIT is a greedy incremental algorithm that starts with the longest input sequence as the first cluster representative and then processes the remaining sequences from long to short to classify each sequence as a redundant or representative sequence based on its similarities to the existing representatives (Fu et al., 2012).


```
cd-hit -i pdb_kunitz_customreported.fasta -o pdb_kunitz_customreported.clstr -c 0.9

```
The files were analyse to check if some clusters needed to be cut away, before proceeding with  Multiple Sequence Alignment. results file: pdb_kunitz_customreported_trimmed.fasta. The most representative ID from each cluster were extracted. 

```

clstr2txt.pl pdb_kunitz_customreported_filtered.clstr > pdb_kunitz.clusters.txt

```
The ids from that file are taken and stored in a FASTA file, and prepare a file for the MSA

```

awk '$5 == 1 {print $1}' pdb_kunitz.clusters.txt > pdb_kunitz_rp.ids

for i in $(cat pdb_kunitz_rp.ids); do   grep -A 1 "^>$i" pdb_kunitz_customreported.fasta | head -n 2 >> pdb_kunitz_rp.fasta; done

```

pdb_kunitz_rp.fasta contains only the represantative sequences selected from CD-HIT

### 2. *Perform MSA*

PDBeFold Multi alignment Tool (https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver) was used. 

after the alignment the file was saved as .ali and then formatted

```
awk '{if (substr($1,1,1)==">") {print "\n" toupper($1)} else {printf "%s", toupper($1)}}'pdb_kunitz_rp.ali >pdb_kunitz_rp_formatted.ali
```

### 3. *Bulid the HMM*

The dedicated conda enviroment was activated and the HMM model was constructed using the output file from our previous MSA. 

```
conda activate 
hmmbuild structural_model.hmm pdb_kunitz_rp_formatted.ali

```

### 4. *Model testing*

to test the model a 2-k fold cross validation was used, creating two positive and negative set. 
The positive set was created removing from a file containing all kunitz proteins ids (downloaded from Uniprot) of the proteins we use to create the model, eliminating any possible bias. 
From uniprot three files were downloaded:
human_kunitz.fasta
human_not_kunitz.fasta
uniprot_sprot.fasta
and then a BLAST search was run.  
```
cat human_kunitz.fasta human_not_kunitz.fasta > all_kunitz.fasta

makeblastdb -in all_kunitz_uniprot.fasta -dbtype prot -out all_kunitz_uniprot.fasta

blastp -query pdb_kunitz_rp.fasta -db all_kunitz_uniprot.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
```
The ids of all the proteins with high identity to the one with witch we have created our HMM model were retrived. 

```
grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > to_remove.ids

grep ">" all_kunitz_uniprot.fasta | cut -d "|" -f 2 > all_kunitz.id

comm -23 <(sort all_kunitz.id) <(sort to_remove.ids) >to_keep.ids

```

### 5. *Format results*

wThe script *get_seq.py* was used to filter the files and prepare positive and negative set. 

```
python3 get_seq.py to_keep.ids all_kunitz_uniprot.fasta ok_kunitz.fasta

grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 >sp.id

comm -23 <(sort sp.id) <(sort all_kunitz.id) >sp_negs.ids

python3 get_seq.py sp_negs.ids uniprot_sprot.fasta sp_negs.fasta

```
The positive and negative set were randomised.

```
sort -R sp_negs.ids > random_sp_negs.ids
sort -R to_keep.ids > random_ok_kunitz.ids

head -n 184 random_ok_kunitz.ids > pos_1.ids
tail -n 184 random_ok_kunitz.ids > pos_2.ids

head -n 286417 random_sp_negs.ids >neg_1.ids
tail -n 286417 random_sp_negs.ids >neg_2.ids

```
Sequences were then retrived using the entire uniprot database.

```
python3 get_seq.py pos_1.ids uniprot_sprot.fasta pos_1.fasta
python3 get_seq.py pos_2.ids uniprot_sprot.fasta pos_2.fasta
python3 get_seq.py neg_1.ids uniprot_sprot.fasta neg_1.fasta
python3 get_seq.py neg_2.ids uniprot_sprot.fasta neg_2.fasta
```

### 6. *Model Evaluation*

Hmmsearch was used, the input was the model and our positives/negatives sets. 

```
hmmsearch -Z 1000 --max --tblout pos_1.out structural_model.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out structural_model.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out structural_model.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out structural_model.hmm neg_2.fasta

```
The output were formatted in a classification format.

```

grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_2.class

```
Two different sets were created from the first four files. 

```
cat pos_2.class neg_2_hits.class >set_2.class
cat pos_1.class neg_1_hits.class > set_1.class
```
The file performance.py was used to value the performance.

```

python3 performance.py set_2.class 1e-5 >results_set_1.txt
python3 performance.py set_2.class 1e-5 >results_set_2.txt

```

Performance.py was run using different thresholds for the E-value, analysing the changes in the performance metrics. 

```

for i in $(seq 1 10); do   python3 performance.py set_1.class 1e-$i; done | sort -nrk 6 > diff_threshold_set1.txt
for i in $(seq 1 10); do   python3 performance.py set_2.class 1e-$i; done | sort -nrk 6 > diff_threshold_set2.txt

```

### 7. Output

 All final evaluation metrics, confusion matrices, and plots are saved in the results/ folder.

### 8. Author

Cagnini Luca
MSc Student in Bioinformatics,
University of Bologna

### License

This project is licensed under the MIT License. See the LICENSE file for details.


