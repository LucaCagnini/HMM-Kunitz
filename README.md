
# Profile HMM Pipeline for Kunitz-Type Protease Inhibitor Domain (Pfam: PF00014)

This repository implements a **computational pipeline** for building and evaluating a Profile Hidden Markov Model (HMM) designed to detect the **Kunitz-type protease inhibitor domain** (Pfam ID: PF00014) within protein sequences.

The project was developed as part of the _Laboratory of Bioinformatics 1_ course (MSc in Bioinformatics, University of Bologna) and integrates **structural bioinformatics**, **sequence analysis**, and **statistical model evaluation**.

---
## Table of Contents

- [Objectives](#objectives)
- [Requirements](#requirements)
- [Pipeline Execution](#pipeline-execution)
  - [1. Extract Kunitz Sequences](#1-extract-kunitz-domain-sequences-from-uniprot)
  - [2. Build the HMM](#2-build-the-hmm)
  - [3. Search with the Model](#3-search-with-the-model)
  - [4. Format Results](#4-format-results-into-class-files)
  - [5. Evaluate the Model](#5-evaluate-the-model)
- [Output](#-output)
- [Author](#-author)
- [License](#-license)

  ---

## objectives

- Build a structure-based multiple sequence alignment of known Kunitz-domain proteins.
- Create and calibrate an HMM model using that alignment.
- Evaluate the model’s performance using curated positive (true Kunitz) and negative (non-Kunitz) datasets.
- Filter redundancy using BLAST and CD-HIT.
- Output metrics such as confusion matrix, accuracy, sensitivity, and specificity.

---

## Requirements

all the requirements needs to be satisfied and can be installed in a dedicated conda enviroment.

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
From the advance search function in UniProt it is possible to download the sequences of all the proteins showing a Kunitz domain:
   
   ° Pfam ID: PF 00014
   
   ° Resolution ≤ 3.5 Angstrom
   
   ° 45 ≤ length ≥ 80
   
Query:

```
Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
```

Then we can extract the IDs from the results file (rcsb_pdb_custom_report20250410062557.csv) :

```
cat rcsb_pdb_custom_report20250410062557.csv| tr -d '"' | awk -F ',' '{if (length($2)>0) {name=$2}; print name ,$3,$4,$5}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz_customreported.fasta
```

Then, we cluster the sequences using CD-HIT, at 90% identity threshold: 

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
 1543  cat pdb_kunitz_rp.ids

```

and prepare a file for the MSA

```
for i in $(cat pdb_kunitz_rp.ids); do   grep -A 1 "^>$i" pdb_kunitz_customreported.fasta | head -n 2 >> pdb_kunitz_rp.fasta; done

```

pdb_kunitz_rp.fasta contains only the represantative sequences selected from CD-HIT

### Perform MSA

We used the PDBeFold Multi alignment Tool
[link text] (https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver)

we download the PDEBefold alignment as .ali and then format that 

```
awk '{if (substr($1,1,1)==">") {print "\n" toupper($1)} else {printf "%s", toupper($1)}}'pdb_kunitz_rp.ali >pdb_kunitz_rp_formatted.ali
```

### Bulid the HMM

We activate our dedicated conda enviroment and build the HMM model using the output file from the previous MSA

```
conda activate 
hmmbuild structural_model.hmm pdb_kunitz_rp_formatted.ali

```

### Search with the model 
We run hmmsearch on all positive and negative FASTA files and create a tabular output. We generate .out files that contains the E-values computed on the HMM of the sequence alignment.

```
hmmsearch -Z 1000 --max --tblout pos_1_seqali.out pdb_kunitz_rp_seqali.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2_seqali.out pdb_kunitz_rp_seqali.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1_seqali.out pdb_kunitz_rp_seqali.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2_seqali.out pdb_kunitz_rp_seqali.hmm neg_2.fasta
```





