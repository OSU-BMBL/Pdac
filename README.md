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









---

## Typical Installation Time
Steps 1–3 (conda environment, PyTorch, PyG extensions and EmitGCL) take about **2 minutes** with the GPU build (`torch==1.12.0+cu113`, prebuilt PyG wheels), and about **15 minutes** with the CPU build, where `torch_scatter`, `torch_sparse` and `torch_cluster` are compiled from source. Download of the input `.h5ad` (about 9.9 GB) is not included.

Reference environment: BMBL server (not OSC), AMD EPYC 7H12 (128 cores), ~1 TB RAM, NVIDIA A100 80GB, RHEL 8.10, Python 3.8, torch 1.12.0+cu113, torch_geometric 2.1.0, EmitGCL 0.2.1.

> Note: the GPU command in step 2 uses `cu102`, which does not support A100 GPUs. We used `cu113` (`pip install torch==1.12.0+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html` and `-f https://data.pyg.org/whl/torch-1.12.0%2Bcu113.html` for the PyG extensions).

## Expected Demo Runtime
Demo command: `python Pancreas-pn0325.py` with `scRNA_harmony.h5ad` in the working directory (step 5).

- **Dataset:** `scRNA_harmony.h5ad`, 77,240 cells × 31,915 genes; 45,312 cells (primary and metastatic samples) are used for training.
- **Workflow:** loading and batch partitioning, NodeDimensionReduction training (10 epochs), EmitGCL training (10 epochs), then prediction and attention output.
- **Runtime:** about **20–25 hours** on one A100 GPU. Loading and batch partitioning take about 20 minutes, and each training epoch takes about 1 hour (measured for the first epoch). The total is this per-epoch time over the 20 epochs.
- **Resources:** 1 GPU (~2.6 GB GPU memory), ~11 GB RAM, and multiple CPU cores (batch partitioning uses 48 worker processes). Training is largely CPU-bound, so the GPU gives little speed-up over the CPU build.

These times were obtained on the BMBL server described above. They have not been measured on OSC or on a desktop computer.
