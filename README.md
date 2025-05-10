
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






