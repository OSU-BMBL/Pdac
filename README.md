# Senescence–EMT Axis Underlies Pancreatic Cancer Progress

## Overview
Pancreatic ductal adenocarcinoma (PDAC) is a highly lethal malignancy driven by early metastatic dissemination, profound immune evasion, and poor therapeutic responsiveness. Existing molecular classifications, largely based on bulk transcriptomic data, have limited capacity to resolve the functional heterogeneity of tumor cells or to reveal actionable mechanisms of progression.  

In this project, we establish a **deep learning–based framework** integrating large-scale single-cell RNA-seq to uncover a **senescence–EMT axis** that stratifies PDAC cells into four discrete tumor states. These states capture distinct transcriptional programs, immune microenvironments, and metastatic potentials, spanning a trajectory from immune-cold (SENE-EMT-) to immune-evasive (SENE-EMT-) phenotypes.

---

## Key Findings
- **New classification of PDAC tumor states**: four discrete states along the senescence–EMT axis.  
- **Clinical significance**: patients enriched for the SENE-EMT- state exhibited significantly better prognosis than all other states.  
- **Cross-platform validation**: confirmed across Visium spatial transcriptomics, CosMX high-plex spatial imaging, and an independent clinical cohort of **34,138 patients** with treatment and follow-up data.  
- **Mechanistic insights**: identified state-specific transcriptional regulators.  
- **Therapeutic mapping**: linked regulators to **FDA-approved drugs** (e.g., metformin, statins, ATRA), providing immediate translational opportunities.  

---

## Graphical Abstract
*<img width="2825" height="1643" alt="image" src="https://github.com/user-attachments/assets/da259cf5-afb3-42cb-8a30-f91901d114ed" />*  

---

## Repository Contents
- `notebooks/` – Jupyter or R notebooks for data preprocessing and analysis.  
- `src/` – deep learning framework implementation.  
- `data/` – metadata or links to raw data (scRNA-seq, Visium, CosMX).  
- `figures/` – main and supplementary figures.  

---

## Installation
Clone this repository and set up the environment:

```bash
git clone https://github.com/YourUsername/Senescence-EMT-PDAC.git
cd Senescence-EMT-PDAC

# Create environment (example with conda)
conda create -n pdac python=3.10
conda activate pdac
pip install -r requirements.txt

