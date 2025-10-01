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
The installation process involves some optional and necessary steps. Here's the detailed breakdown:

### 1) Create a new environment
```bash
conda create -n pdac python=3.8
conda activate pdac
```

### 2）PyTorch install
#### CPU version (For Linux/Windows/MacOS system (torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):)
```bash
pip install torch==1.12.0+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html
pip install torch_scatter==2.0.9 torch_sparse==0.6.14 torch_cluster==1.6.0 -f https://data.pyg.org/whl/torch-1.12.0%2Bcpu/
```
#### GPU version
Please visit the official PyTorch website at PyTorch to select and download the CUDA-enabled version of PyTorch that best matches your system configuration.
For linux system(You need to select the version that is compatible with your system's graphics card. For example: torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.1.0+ torch_sparse-0.6.16):
```batch
pip install torch==1.12.0+cu102 -f https://download.pytorch.org/whl/cu102/torch_stable.html
pip install torch_scatter==2.1.0 torch_sparse==0.6.16 torch_cluster==1.6.0 -f https://data.pyg.org/whl/torch-1.12.0%2Bcu102/
```
### 3）Install EmitGCL
```batch
pip install --upgrade EmitGCL
```
### 4）Data enhancement
Prepare matched scRNA-seq data from the primary and metastatic site.
Prepare scRNA-seq data in one of two sites.
Run the Pancreas-matched.R to get the processed .h5ad.

### 5）Model running
Input: .h5ad file from the data enhancement step. You can also download it from [here](https://drive.google.com/drive/folders/13enekPdeDxLEAkGK4rrXkeBkGtlKX3ZT?usp=sharing).

```python
cd Pdac & Pancreas-pn0325.py
```

Output: cell embedding, predicted cell clusters.








