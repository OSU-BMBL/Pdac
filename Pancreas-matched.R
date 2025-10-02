setwd('/yourPath/Pdac/')
# Load necessary libraries
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(ggplot2)
library(reticulate)
library(harmony)
library(dplyr)
library(patchwork)

# Import Python 'anndata' package
ad <- import("anndata")

# Set the patient directories
patients <- c("Patient1", "Patient2", "Patient3")
base_dir <- "/yourPath/GSE197177"

# Initialize list to store Seurat objects for each patient
seurat_list <- list()

# Iterate through each patient's directory and process the data
for (patient in patients) {
  # Set the working directory for the current patient
  patient_dir <- file.path(base_dir, patient)
  setwd(patient_dir)
  
  # Define file paths for data files
  cell_label_file <- "cell_label.csv"       # Cell label file
  cell_barcodes_file <- "cell_barcodes.csv" # Cell barcode file
  gene_cell_file <- "gene_cell.mtx"         # Gene expression matrix file
  gene_symbol_file <- "gene_symbol.csv"     # Gene symbol file
  
  # Read the data files
  cell_label <- read.csv(cell_label_file, header = TRUE, stringsAsFactors = FALSE)  
  cell_barcodes <- read.csv(cell_barcodes_file, header = FALSE, stringsAsFactors = FALSE)[-1, ] 
  gene_cell <- readMM(gene_cell_file)  
  gene_symbol <- read.csv(gene_symbol_file, header = FALSE, stringsAsFactors = FALSE)[-1, ]  
  
  # Check for consistency between cell_label and cell_barcodes
  if (nrow(cell_label) != length(cell_barcodes)) {
    stop(paste("cell_label and cell_barcodes row count mismatch for patient", patient))
  }
  
  # Create unique cell names by combining Pid and barcode
  cell_names <- paste(cell_label$Pid, cell_barcodes, sep = "_")
  
  # Ensure unique gene names and set cell names as column names
  rownames(gene_cell) <- make.unique(gene_symbol)
  colnames(gene_cell) <- cell_names
  rownames(cell_label) <- cell_names
  
  # Add patient information to the metadata
  cell_label$Patient <- patient  
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = gene_cell, meta.data = cell_label)
  
  # Add Seurat object to the list
  seurat_list[[patient]] <- seurat_obj
}

# Merge all patients' Seurat objects into one combined object
combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = patients)

# Check if "Pid" column exists and filter out cells with "Pid" = "N"
if ("Pid" %in% colnames(combined_seurat@meta.data)) {
  combined_seurat <- subset(combined_seurat, subset = Pid != "N")
  cat("Number of cells after filtering Pid = 'N':\n")
  print(dim(combined_seurat))
} else {
  cat("$Pid column is missing, please check metadata!\n")
}

# Create a new column 'Patient_Pid' by merging 'Patient' and 'Pid'
combined_seurat$Patient_Pid <- paste(combined_seurat$Patient, combined_seurat$Pid, sep = "_")

# Filter files containing 'PN' in the path
file_path <- "/yourPath"
rds_file_paths <- c()

if (dir.exists(file_path)) {
  files_and_dirs <- list.files(file_path, full.names = TRUE)
  filtered_files <- files_and_dirs[!grepl("2724|2760", basename(files_and_dirs))]
  for (subdir in filtered_files) {
    if (dir.exists(subdir)) {
      rds_in_subdir <- list.files(subdir, pattern = "\\.rds$", full.names = TRUE)
      if (length(rds_in_subdir) > 0) {
        rds_file_paths <- c(rds_file_paths, rds_in_subdir[grepl("PN", rds_in_subdir)])
      }
    }
  }
  cat("Filtered .rds file paths containing 'PN':\n")
  print(rds_file_paths)
} else {
  cat("The directory does not exist, please check the path: ", file_path, "\n")
}

# Iterate over filtered .rds files and process each
for (rds_file_path in rds_file_paths) {
  # Extract the second-to-last part of the path as the directory name
  dir_name <- strsplit(rds_file_path, "/")[[1]]
  dir_name <- dir_name[length(dir_name) - 1]
  
  # Set the save path
  save_path <- file.path("/yourPath", dir_name)
  
  # Create the save path if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Check if the file exists and process it
  if (file.exists(rds_file_path)) {
    # Load Seurat object from the .rds file
    seurat_obj <- readRDS(rds_file_path)
    cat("Processing file:", rds_file_path, "\n")
    
    # Get common genes between the current Seurat object and combined Seurat object
    genes_seurat_obj <- rownames(seurat_obj@assays$RNA)
    genes_combined <- rownames(combined_seurat@assays$RNA)
    common_genes <- intersect(genes_seurat_obj, genes_combined)
    cat("Number of common genes:", length(common_genes), "\n")
    
    # Subset both Seurat objects to only include the common genes
    seurat_obj <- subset(seurat_obj, features = common_genes)
    combined_seurat <- subset(combined_seurat, features = common_genes)
    
    # Add Patient_Pid column to the Seurat object
    seurat_obj$Patient_Pid <- paste0(dir_name, "_P")
    
    # Merge the current Seurat object with the combined Seurat object
    merged_seurat <- merge(combined_seurat, y = seurat_obj, add.cell.ids = c("Combined", dir_name))
    
    # Perform Harmony batch correction
    scRNA_harmony <- merged_seurat %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(verbose = FALSE)
    
    scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "Patient_Pid")
    scRNA_harmony <- scRNA_harmony %>%
      RunUMAP(reduction = "harmony", dims = 1:30) %>%
      FindNeighbors(reduction = "harmony", dims = 1:30) %>%
      FindClusters()
    
    # Plot UMAP and save the image
    plot1 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
    plot2 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "Patient_Pid")
    combined_plot <- plot1 | plot2
    output_file <- file.path(save_path, "combined_umap_plot.png")
    ggsave(filename = output_file, plot = combined_plot, width = 12, height = 6)
    cat("UMAP plot saved to:", output_file, "\n")
    
    # Save Seurat object as h5Seurat format
    h5seurat_file <- file.path(save_path, "scRNA_harmony.h5Seurat")
    SaveH5Seurat(
      scRNA_harmony,
      filename = h5seurat_file,
      overwrite = TRUE,
      assay = "RNA",
      feature.data = FALSE,
      scale.data = FALSE,
      raw.data = TRUE,
      features = "all"
    )
    
    # Load h5Seurat file and convert to h5ad format
    scRNA_harmony_loaded <- LoadH5Seurat(h5seurat_file)
    counts_matrix <- t(as.matrix(scRNA_harmony_loaded@assays$RNA@counts))
    obs_data <- scRNA_harmony_loaded@meta.data
    var_data <- data.frame(genes = rownames(scRNA_harmony_loaded))
    
    adata <- ad$AnnData(X = counts_matrix, obs = obs_data, var = var_data, dtype = "float32")
    h5ad_file <- file.path(save_path, "scRNA_harmony.h5ad")
    adata$write_h5ad(h5ad_file)
    cat("Seurat object saved as h5ad format, file path:", h5ad_file, "\n")
    
  } else {
    cat("File does not exist, please check the path:", rds_file_path, "\n")
  }
}

####
sessionInfo()
#############
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 9.4 (Plow)
#
# Matrix products: default
# BLAS/LAPACK: /apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-mkl/gcc/11.4.1/2023.2.0-gwnin2p/mkl/2023.2.0/lib/intel64/libmkl_gf_lp64.so.2;  LAPACK version 3.10.1
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#
# time zone: US/Eastern
# tzcode source: internal
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
# [1] future_1.67.0         patchwork_1.3.2       dplyr_1.1.4           harmony_1.2.3        
# [5] Rcpp_1.1.0            reticulate_1.43.0     ggplot2_3.5.2         Matrix_1.7-4         
# [9] SeuratDisk_0.0.0.9021 Seurat_5.3.0          SeuratObject_5.2.0    sp_2.2-0             

