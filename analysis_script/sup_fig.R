library(CellChat)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs)
library(tidyr)
library(purrr)
library(stringr)
rm(list = ls())
setwd('/bmbl_data/weidong/xy/sup/')
set.seed(1111)
sene_marker <- c('LMNB1','CCL4','CCL5','CXCL1','CXCL10','IL1A','IL6','TNF','MMP2','SERPINE1',
                 'IGFBP3','ANKRD1','CCL11','CD68','CXCL12','ICAM1','HMGB1','SMAD4','TP53','KRAS')
emt_marker  <- c("ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2","VIM","FN1","CDH2","ITGB6",
                 "MMP2","MMP9","MMP14","SPARC","COL1A1","COL3A1","COL5A1","TNC","PDGFRB",
                 "TGFBI","FAP","LAMC2","LGALS1")
cancer_marker <- c('EPCAM','KRT19','KRT7','CD44')

classical_genes <- c(
  "GATA6", "PDX1", "HNF1A", "HNF4A", "FOXA2", "TFF1", "MUC1",
  "AGR2", "SPDEF", "EPCAM", "ERBB3", "CDH1", "CA2", "CLDN4"
)
basal_genes <- c(
  "TP63", "KRT5", "KRT14", "S100A2", "ZEB1", "SNAI2", "VIM",
  "CXCL1", "JUN", "MYC", "CDH2", "KRT17", "FOSL1"
) 



folder_path <- "/bmbl_data/xy/plotqs"
# List all .qs files in that folder

cancer <- readRDS("/bmbl_data/xy/plotqs/EMT.rds")
cancer[[1]]

qs_files <- list.files(folder_path, pattern = "\\.qs$", full.names = TRUE)

cluster = qread('/bmbl_data/xy/plotqs/plots_cluster.qs')
cluster[['IU_PDA_HM9']]

# Load them all into a named list
plots_data <- lapply(qs_files, qread)

plots_data[[1]]




####################below ###################
# 肿瘤/免疫/基质的大类注释
auto_pt_size <- function(obj, image_name) {
  sf <- obj@images[[image_name]]@scale.factors
  if (!is.null(sf$fiducial) && !is.null(sf$spot)) {
    return(sf$fiducial / sf$spot)
  } else {
    return(1)
  }
}
markers_list <- list(
  Tumor      = c("EPCAM","KRT8","KRT18","KRT19","MSLN","MUC1","CLDN4","CEACAM6","TSPAN8","AGR2"),
  CAF        = c("COL1A1","COL1A2","DCN","LUM","FAP","PDPN","TAGLN","ACTA2","COL3A1","SPARC","MMP2"),
  Endothelial= c("PECAM1","VWF","KDR","PLVAP","EMCN"),
  Myeloid    = c("LYZ","C1QB","C1QA","C1QC","CD68","ITGAM","MSR1","MRC1","CSF1R","S100A8","S100A9"),
  Tcell      = c("CD3D","CD3E","CD2","CD247","CD8A","CD8B","PDCD1","LAG3","TIGIT"),
  Bcell      = c("MS4A1","CD79A","CD79B","MZB1","JCHAIN"),
  NK         = c("NKG7","GNLY","PRF1","GZMB"),
  cDC1       = c("CLEC9A","XCR1","BATF3","IRF8","HLA-DRA","HLA-DRB1"),
  cDC2       = c("CD1C","CLEC10A","FCER1A","ITGAX","HLA-DPB1","HLA-DPA1")
)

SASP_genes <- c("IL6","CXCL8","CXCL1","CXCL2","CCL20","SERPINE1","ICAM1","MMP1","MMP3","TNFRSF10B")
TGFb_targets <- c("TGFB1","TGFBR1","TGFBR2","SMAD3","SMAD4","SERPINE1","CTGF","TAGLN","COL1A1","COL1A2")
TLS_genes  <- c("CXCL13","CCL19","CCL21","LTB","MS4A1","CD79A","ICOS","CXCR5")   # TLS 代表性集合
cDC1_genes <- markers_list$cDC1
library(Seurat)
library(patchwork)


# 1) 选择要画的 IMAGE 名（如果你已知道就直接写 IMAGE <- "你的图名"）

######
# Define a list of cluster mappings
cluster_map <- list(
  IU_PDA_HM13 = c ("3"),
  IU_PDA_HM2   = c("4","2","3"),
  IU_PDA_HM2_2 = c("2"),
  IU_PDA_T2    = c("8","3"),
  IU_PDA_T1    = c("6"),        # or "6" #T1 1.31 (if you want both, you can store as list)
  IU_PDA_T3    = c("7","3","12"),
  IU_PDA_HM11  = c("10","11","4","0","12","6"),
  IU_PDA_T11   = c("4","10"),
  IU_PDA_HM9   = c("7"),
  IU_PDA_T9    = c("7"),
  IU_PDA_LNM7  = c("10"),
  IU_PDA_HM4   = c("3","4"),
  IU_PDA_T4    = c("3","10","8","5","2","4"),
  IU_PDA_HM12  = c("10","5","4","1","3","9","0","7","2","6"),
  IU_PDA_T12   = c("3","5","4","11","9","0","7","2"),
  IU_PDA_LNM12 = c("3","4","6","8"),
  IU_PDA_HM6   = c("4"),
  IU_PDA_LNM6  = c("9","4","7","10","2","0","6","11"),
  IU_PDA_T6    = c("5","4","8","3","2","10"),
  IU_PDA_HM10  = c("2"),
  IU_PDA_LNM10 = c("0","1","2","3","4","5","6","7","8","9","10"),
  IU_PDA_T10   = c("9","3"),
  IU_PDA_T8    = c("3","8","9","5"),
  IU_PDA_LNM8  = c("1","4","9","5","2"),
  IU_PDA_HM8   = c("4")
)

# Define a function to retrieve clusters by sample name
get_clusters <- function(sample_name) {
  if (sample_name %in% names(cluster_map)) {
    return(cluster_map[[sample_name]])
  } else {
    warning(paste("Sample", sample_name, "not found in cluster_map"))
    return(NULL)
  }
}

# Example usage:

# [1] "6"
imgs <- names(brain@images)

# Extract person ID: digits immediately following HM/T/LNM/NP
get_person_id <- function(x) {
  pat <- "^.*_(HM|T|LNM|NP)([0-9]+).*"
  if (grepl(pat, x)) sub(pat, "\\2", x) else NA_character_
}

person_id <- vapply(imgs, get_person_id, character(1))
img_groups <- split(imgs, person_id)

for (pid in names(img_groups)[6:length(img_groups)]) {
  pdf(sprintf("/bmbl_data/weidong/xy/sup/ExtendFig5.figure/Person_%s_plots.pdf", pid), width = 14, height = 6, onefile = TRUE)
  for (IMAGE in img_groups[[pid]]) {
    print(IMAGE)
    if (grepl("NP", IMAGE)) next 
    message(">> Processing image: ", IMAGE)
    
    cells_to_keep <- rownames(brain@images[[IMAGE]]@coordinates)
    counts_mat <- brain@assays$Spatial@counts
    colnames(counts_mat) <- rownames(brain@meta.data)          
    counts_mat_subset <- counts_mat[, cells_to_keep, drop = FALSE]
    meta_data_subset  <- brain@meta.data[cells_to_keep, , drop = FALSE]
    image_subset      <- brain@images[[IMAGE]]
    rownames(image_subset@coordinates) <- rownames(meta_data_subset)
    
    brain_subset <- CreateSeuratObject(counts = counts_mat_subset,
                                       meta.data = meta_data_subset,
                                       assay = "Spatial")
    brain_subset@images[[IMAGE]] <- image_subset
    
    pt_factor <- auto_pt_size(brain_subset, IMAGE)*1.2
    DefaultAssay(brain_subset) <- "Spatial"
    brain_subset <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
    DefaultAssay(brain_subset) <- "SCT"
    brain_subset <- RunPCA(brain_subset, verbose = FALSE)
    brain_subset <- FindNeighbors(brain_subset, reduction = "pca", dims = 1:30)
    brain_subset <- FindClusters(brain_subset, graph.name = "SCT_snn", resolution = 0.6, verbose = FALSE)
    brain_subset <- RunUMAP(brain_subset, reduction = "pca", dims = 1:30)
    
    
    
    all_clusters <- levels(brain_subset$seurat_clusters)
    cluster_colors <- setNames(rep("lightgray", length(all_clusters)), all_clusters)
    cancer_clusters = get_clusters(IMAGE)
    cluster_colors[cancer_clusters] <- "#F5A5BE"
    
    vln_by_cluster <- function(obj, features, group.by = "seurat_clusters", assay = "SCT",
                               ncol = 4, same.y.lims = TRUE, pt.size = 0.0) {
      feats <- intersect(features, rownames(obj))
      if (length(feats) == 0) stop("no in obj：", paste(features, collapse=", "))
      if (length(feats) < length(features)) {
        message("skiped genes： ",
                paste(setdiff(features, feats), collapse=", "))
      }
      plots <- VlnPlot(
        obj, features = feats, group.by = group.by, assay = assay,
        pt.size = pt.size, combine = FALSE, same.y.lims = same.y.lims
      )
      plots <- lapply(plots, function(p)
        p + theme_classic(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                legend.position = "none") +
          xlab(NULL) + ylab("log1p expr"))
      wrap_plots(plots, ncol = ncol)
    }
    p_myeloid <- vln_by_cluster(brain_subset, markers_list$Myeloid, group.by = "seurat_clusters", ncol = 4)
    p_myeloid
    
    pt_factor <- auto_pt_size(brain_subset, IMAGE)*1.2
    SpatialDimPlot(brain_subset,group.by="seurat_clusters",)
    
    
    pt_factor <- auto_pt_size(brain_subset, IMAGE)*1.1
    c <- SpatialDimPlot(
      brain_subset,
      group.by = "seurat_clusters",
      pt.size.factor = pt_factor,
      images = IMAGE,
    cols = cluster_colors
    )
    

    CAF_barrier   <- c("FAP","PDPN","TAGLN","ACTA2","COL1A1","COL1A2","SPARC","DCN","LUM")
    SASP_genes    <- c("IL6","CXCL8","CCL20","MMP1","MMP3","SERPINE1","ICAM1")
    TGFb_targets  <- c("TGFB1","TGFBR1","TGFBR2","SMAD3","SMAD4","SERPINE1","CTGF")
    cDC1_set      <- c("CLEC9A","XCR1","BATF3","IRF8","HLA-DRA","HLA-DRB1")
    TLS_set       <- c("CXCL13","CCL19","CCL21","LTB","MS4A1","CD79A","BCL6","ICOS")
    Treg_set      <- c("FOXP3","IL2RA","CTLA4","TNFRSF18")
    MDSC_set      <- c("S100A8","S100A9","LST1","OLR1","VCAN","CCR2","LILRB4")
    Exhaust_set   <- c("PDCD1","LAG3","TIGIT","HAVCR2")
    CXCL12_axis   <- c("CXCL12","CXCR4")
    CCL2_axis     <- c("CCL2","CCR2")
    CSF1_axis     <- c("CSF1","CSF1R")
    
    #add_score <- function(obj, gs, nm) AddModuleScore(obj, features=list(intersect(gs, rownames(obj))), name=paste0(nm,"_"))
    add_score <- function(obj, gs, nm, nbin = 24) {
      tryCatch(
        {
          AddModuleScore(
            obj,
            features = list(intersect(gs, rownames(obj))),
            name = paste0(nm, "_"),
            nbin = nbin
          )
        },
        error = function(e) {
          message("⚠️ Error with nbin = ", nbin, " → retrying with nbin = 12")
          AddModuleScore(
            obj,
            features = list(intersect(gs, rownames(obj))),
            name = paste0(nm, "_"),
            nbin = 12
          )
        }
      )
    }
    
    
    brain_subset <- add_score(brain_subset, CAF_barrier, "CAF",nbin=12)
    brain_subset <- add_score(brain_subset, SASP_genes, "SASP")
    brain_subset <- add_score(brain_subset, TGFb_targets, "TGFb")
    brain_subset <- add_score(brain_subset, cDC1_set, "cDC1")
    brain_subset <- add_score(brain_subset, TLS_set,  "TLS")
    brain_subset <- add_score(brain_subset, Treg_set, "Treg")
    brain_subset <- add_score(brain_subset, MDSC_set,"MDSC")
    brain_subset <- add_score(brain_subset, Exhaust_set,"Exh")
    brain_subset <- add_score(brain_subset, CXCL12_axis,"CXCL12ax")
    brain_subset <- add_score(brain_subset, CCL2_axis,  "CCL2ax")
    brain_subset <- add_score(brain_subset, CSF1_axis,  "CSF1ax")
    
    
    d <- SpatialFeaturePlot(brain_subset, features=c("CAF_1","CXCL12ax_1","SASP_1","TGFb_1","cDC1_1","TLS_1"), 
                            ncol = 6,pt.size.factor=pt_factor)+ 
      plot_annotation(title = IMAGE)
    
    print(d)
  }
  dev.off()
}


####################above ###################

############immune_Supr ##################


immune_suppression_high <- c(
  "CD274", "CD47", "LGALS9", "IDO1", "IDO2", "VTCN1", "CD276", "FASLG", "TNFSF10",
  "TGFB1", "IL10", "VEGFA", "CCL2", "CXCL8", "AXL", "TYRO3", "CTNNB1", "C10orf54", "SERPINB9", "LGALS1",
  "SOCS1", "SOCS3"
)
for (pid in names(img_groups)[1:length(img_groups)]) {
  pdf(sprintf("/bmbl_data/weidong/xy/sup/immune_suppression/Person_%s_plots.pdf", pid), width = 6, height = 6, onefile = TRUE)
  for (IMAGE in img_groups[[pid]]) {
    print(IMAGE)
    if (grepl("NP", IMAGE)) next 
    message(">> Processing image: ", IMAGE)
    
    cells_to_keep <- rownames(brain@images[[IMAGE]]@coordinates)
    counts_mat <- brain@assays$Spatial@counts
    colnames(counts_mat) <- rownames(brain@meta.data)          
    counts_mat_subset <- counts_mat[, cells_to_keep, drop = FALSE]
    meta_data_subset  <- brain@meta.data[cells_to_keep, , drop = FALSE]
    image_subset      <- brain@images[[IMAGE]]
    rownames(image_subset@coordinates) <- rownames(meta_data_subset)
    
    brain_subset <- CreateSeuratObject(counts = counts_mat_subset,
                                       meta.data = meta_data_subset,
                                       assay = "Spatial")
    brain_subset@images[[IMAGE]] <- image_subset
    
    pt_factor <- auto_pt_size(brain_subset, IMAGE)*1.2
    DefaultAssay(brain_subset) <- "Spatial"
    brain_subset <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
    DefaultAssay(brain_subset) <- "SCT"
    brain_subset <- RunPCA(brain_subset, verbose = FALSE)
    brain_subset <- FindNeighbors(brain_subset, reduction = "pca", dims = 1:30)
    brain_subset <- FindClusters(brain_subset, graph.name = "SCT_snn", resolution = 0.6, verbose = FALSE)
    brain_subset <- RunUMAP(brain_subset, reduction = "pca", dims = 1:30)
    brain_subset$immune_Supr   <- as.numeric(score_by_mean(brain_subset, immune_suppression_high,   assay = "SCT"))
    
    pt_factor <- auto_pt_size(brain_subset, IMAGE)*1.1
    immune_Supr <- SpatialFeaturePlot(
      brain_subset,
      images   = IMAGE,
      features = "immune_Supr",
      pt.size.factor = pt_factor, max.cutoff=0.3
    ) +
      #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
      plot_annotation(title = IMAGE)
    
    
    print(immune_Supr)
  }
  dev.off()
}

############immune_Supr  above ##################

############### net bubble ###############
readRDS('/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp1_full.rds')

netVisual_bubble


############## net bubble above #########
