library(CellChat)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
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


pdac <- readRDS("/bmbl_data/xy/PDAC_Updated.rds")
sample_label<- sub("(^[^_]+_[^_]+).*", "\\1", names(Idents(pdac)))
pdac$sample_label <- sample_label
my_colors <- c(
  "#F2778F", "#F48AAE", "#F6A6C1", "#F8B4C9", "#F39CA8",
  "#FF8E66", "#FFA066", "#FFB26B", "#FFC57A", "#FFD48F",
  "#FFE07A", "#F3E562", "#D6EB6E", "#BDE67A", "#A6E08B",
  "#7FD6A6", "#64CFBA", "#5FD1D1", "#6FD7E3", "#86DFF0",
  "#8CC6F7", "#86A9F6", "#9A95F5", "#B08FF2", "#C88AEF"
)
brain <- pdac

##############################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

## 1) Gene sets (can be modified as needed)
sene_marker <- c('LMNB1','CCL4','CCL5','CXCL1','CXCL10','IL1A','IL6','TNF','MMP2','SERPINE1',
                 'IGFBP3','ANKRD1','CCL11','CD68','CXCL12','ICAM1','HMGB1','SMAD4','TP53','KRAS')
emt_marker  <- c("ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2","VIM","FN1","CDH2","ITGB6",
                 "MMP2","MMP9","MMP14","SPARC","COL1A1","COL3A1","COL5A1","TNC","PDGFRB","TGFBI","FAP","LAMC2","LGALS1")
cancer_marker <- c('EPCAM','KRT19','KRT7','CD44')

## 2) Scoring function (calculate column means from log-normalized matrix with rows=genes, columns=cells)
score_by_mean <- function(seu, genes, assay = DefaultAssay(seu), slot = "data") {
  mat <- GetAssayData(seu, assay = assay, slot = slot)
  if (nrow(mat) == 0) return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  row_up <- toupper(rownames(mat))
  idx <- match(toupper(genes), row_up)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  sc <- Matrix::colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  sc[is.nan(sc)] <- NA_real_
  sc
}

## 3) Colors
cols_state <-c("Sene+EMT+" = "#F0BCBE", 
               "Sene+EMT-" = "#EB7B75", 
               "Sene-EMT+"      = "#F2D2AB",
               "Sene-EMT-"      = "#EEAA7B",
               "Other"     = "grey90" )

genes_of_interest <- c("SOX2", "NANOG", "POU5F1", "KLF4", "MYC")
geneSets <- c("SOX2", "NANOG", "POU5F1", "KLF4", "MYC")
classical_genes <- c(
  "GATA6", "PDX1", "HNF1A", "HNF4A", "FOXA2", "TFF1", "MUC1",
  "AGR2", "SPDEF", "EPCAM", "ERBB3", "CDH1", "CA2", "CLDN4"
)
basal_genes <- c(
  "TP63", "KRT5", "KRT14", "S100A2", "ZEB1", "SNAI2", "VIM",
  "CXCL1", "JUN", "MYC", "CDH2", "KRT17", "FOSL1"
)

immune_suppression_high <- c(
  "CD274", "CD47", "LGALS9", "IDO1", "IDO2", "VTCN1", "CD276", "FASLG", "TNFSF10",
  "TGFB1", "IL10", "VEGFA", "CCL2", "CXCL8", "AXL", "TYRO3", "CTNNB1", "C10orf54", "SERPINB9", "LGALS1",
  "SOCS1", "SOCS3"
)
sample_list <- sample_list[29:46]
geneSets <- list("Stemness_Score" = genes_of_interest)
Fentropy <- function(x) {
  prob <- x / sum(x)
  -sum(prob * log2(prob + 1e-8))
}
plots_state <- list()
plots_cancerHL <- list()
all_cells <- list()
all_comp  <- list()
plots_EMT <- list()
plots_SENE <- list()
plots_stemness <- list()
plots_entropy <- list()
plots_basal <-list()
plots_classical <- list()
plots_cluster <- list()
plots_immune <- list()

IMAGE=names(brain@images)[26]
for (IMAGE in names(brain@images)) {
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
  
  DefaultAssay(brain_subset) <- "Spatial"
  brain_subset <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
  DefaultAssay(brain_subset) <- "SCT"
  brain_subset <- RunPCA(brain_subset, verbose = FALSE)
  brain_subset <- FindNeighbors(brain_subset, reduction = "pca", dims = 1:30)
  brain_subset <- FindClusters(brain_subset, graph.name = "SCT_snn", resolution = 0.6, verbose = FALSE)
  brain_subset <- RunUMAP(brain_subset, reduction = "pca", dims = 1:30)
  brain_subset$classical   <- as.numeric(score_by_mean(brain_subset, classical_genes,   assay = "SCT"))
  brain_subset$basal   <- as.numeric(score_by_mean(brain_subset, basal_genes,   assay = "SCT"))
  brain_subset$cancer <- as.numeric(score_by_mean(brain_subset, cancer_marker, assay = "SCT"))
  brain_subset$immune_Supr   <- as.numeric(score_by_mean(brain_subset, immune_suppression_high,   assay = "SCT"))
  
  lev <- levels(factor(brain_subset$seurat_clusters))  
  pal <- if (length(lev) > length(my_colors)) {
    grDevices::colorRampPalette(my_colors)(length(lev))
  } else {
    my_colors[seq_along(lev)]
  }
  names(pal) <- lev  
  pal <- c(pal,'gray90')
  names(pal) <-c(lev,'Other')
  cl_mean <- brain_subset@meta.data %>%
    as_tibble(rownames = "spot") %>%
    group_by(seurat_clusters) %>%
    summarise(mean_cancer = mean(cancer, na.rm = TRUE), .groups = "drop")
  immune_Supr <- SpatialFeaturePlot(
    brain_subset,
    images   = IMAGE,
    features = "immune_Supr",
    pt.size.factor = 1.6, max.cutoff=0.3
  ) +
    scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  thr <- median(cl_mean$mean_cancer, na.rm = TRUE)  # can also use mean or quantiles
  cancer_clusters <- cl_mean$seurat_clusters[cl_mean$mean_cancer > thr]
  brain_subset$is_cancer <- brain_subset$seurat_clusters %in% cancer_clusters
  this_cells <- brain_subset@meta.data %>%
    as_tibble(rownames = "spot") %>%
    transmute(sample = IMAGE, spot, seurat_clusters,
              is_cancer = is_cancer, cancer = cancer, basal = basal, classical = classical, immune = immune_Supr )
  all_cells[[IMAGE]] <- this_cells
  
  plots_immune[[IMAGE]] <- immune_Supr
}
library(dplyr)
library(tidyr)
library(ggplot2)

# Combine all IMAGE tibbles; keep IMAGE name as a column
df <- bind_rows(all_cells, .id = "IMAGE")

# Filter cancer spots and reshape basal/classical to long format
df_long <- df %>%
  filter(is_cancer == TRUE) %>%
  select(IMAGE, sample, spot, basal, classical) %>%
  pivot_longer(cols = c(basal, classical),
               names_to = "program",
               values_to = "score") %>%
  drop_na(score)

# Boxplot: basal vs classical within each IMAGE
ggplot(df_long, aes(x = program, y = score, fill = program)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Score", title = "Basal vs Classical (cancer spots only)") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  facet_wrap(~ IMAGE, scales = "free_y")   # one panel per IMAGE

ord <-c('IU_PDA_T2','IU_PDA_HM2','IU_PDA_HM2_2','IU_PDA_T3','IU_PDA_HM3',
        'IU_PDA_T4','IU_PDA_HM4','IU_PDA_T6','IU_PDA_LNM6','IU_PDA_HM6',
        'IU_PDA_T8','IU_PDA_LNM8','IU_PDA_HM8', 'IU_PDA_T9','IU_PDA_HM9',
        'IU_PDA_T10','IU_PDA_LNM10','IU_PDA_HM10','IU_PDA_T11','IU_PDA_HM11',
        'IU_PDA_T12','IU_PDA_HM12')
df_long <- df_long[df_long$sample %in% ord,]
df_long$sample <- factor(df_long$sample, levels = ord)

# Plot
ggplot(df_long, aes(x = sample, y = score, fill = program)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3,
               position = position_dodge(width = 0.7)) +
  geom_jitter(aes(color = program),
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.7),
              size = 0.25, alpha = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = c(basal = "#F6A6C1", classical = "#86CFBA"),
                    name = NULL, labels = c("basal", "classical")) +
  scale_color_manual(values = c(basal = "#F6A6C1", classical = "#86CFBA")) +
  labs(x = "Sample", y = "Score",
       title = "Per-sample basal & classical scores (cancer spots)") +
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(
    axis.line   = element_line(color = "black", linewidth = 0.8),  
    axis.ticks  = element_line(linewidth = 1.0),                  
    axis.text   = element_text(size = 13, family = "Arial"),       
    axis.title  = element_text(size = 15, family = "Arial"),        
    plot.title  = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title= element_text(size = 13, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 8,family = "Arial"),
    axis.title.x = element_blank()
  )


IMAGE <- names(brain@images)[17]
for (IMAGE in names(brain@images)) {
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
  

  DefaultAssay(brain_subset) <- "Spatial"
  brain_subset <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
  DefaultAssay(brain_subset) <- "SCT"
  brain_subset <- RunPCA(brain_subset, verbose = FALSE)
  brain_subset <- FindNeighbors(brain_subset, reduction = "pca", dims = 1:30)
  brain_subset <- FindClusters(brain_subset, graph.name = "SCT_snn", resolution = 0.6, verbose = FALSE)
  brain_subset <- RunUMAP(brain_subset, reduction = "pca", dims = 1:30)
  brain_subset$classical   <- as.numeric(score_by_mean(brain_subset, classical_genes,   assay = "SCT"))
  brain_subset$basal   <- as.numeric(score_by_mean(brain_subset, basal_genes,   assay = "SCT"))
  
  brain_subset$sene   <- as.numeric(score_by_mean(brain_subset, sene_marker,   assay = "SCT"))
  brain_subset$immune_Supr   <- as.numeric(score_by_mean(brain_subset, immune_suppression_high,   assay = "SCT"))
  brain_subset$emt    <- as.numeric(score_by_mean(brain_subset, emt_marker,    assay = "SCT"))
  brain_subset$cancer <- as.numeric(score_by_mean(brain_subset, cancer_marker, assay = "SCT"))
  brain_subset$stemness <- as.numeric(score_by_mean(brain_subset, geneSets, assay = "SCT"))
  lev <- levels(factor(brain_subset$seurat_clusters))  
  pal <- if (length(lev) > length(my_colors)) {
    grDevices::colorRampPalette(my_colors)(length(lev))
  } else {
    my_colors[seq_along(lev)]
  }
  names(pal) <- lev  
  pal <- c(pal,'gray90')
  names(pal) <-c(lev,'Other')
  
  
  cl_mean <- brain_subset@meta.data %>%
    as_tibble(rownames = "spot") %>%
    group_by(seurat_clusters) %>%
    summarise(mean_cancer = mean(cancer, na.rm = TRUE), .groups = "drop")
  
  thr <- median(cl_mean$mean_cancer, na.rm = TRUE)*1.25 
  cancer_clusters <- cl_mean$seurat_clusters[cl_mean$mean_cancer > thr]
  #cancer_clusters <- c("4","2","3") #HM2
  #cancer_clusters <- c("2") #HM2_2
  cancer_clusters <- c("8","3") #T2
  brain_subset$is_cancer <- brain_subset$seurat_clusters %in% cancer_clusters
  
  idx <- which(brain_subset$is_cancer)
  s_med <- median(brain_subset$sene[idx], na.rm = TRUE)
  e_med <- median(brain_subset$emt [idx], na.rm = TRUE)
  
  brain_subset$state <- "Other"
  brain_subset$state[idx] <- paste0(
    ifelse(brain_subset$sene[idx] >= s_med, "Sene+", "Sene-"),
    ifelse(brain_subset$emt [idx] >= e_med,  "EMT+",  "EMT-")
  )
  brain_subset$state <- factor(
    brain_subset$state,
    levels = c("Sene+EMT+","Sene+EMT-","Sene-EMT+","Sene-EMT-","Other")
  )
  
  
  #RNA_matrix<-brain_subset@assays$SCT@data
  #RNA_matrix <- t(RNA_matrix)
  #entropy_scores <- apply(RNA_matrix, 1, Fentropy)
  #brain_subset@meta.data$entropy <- entropy_scores
 
  classical_expr <- rowMeans(as.matrix(RNA_matrix[, intersect(classical_genes,rownames(brain_subset))]), na.rm = TRUE)
  basal_expr <- rowMeans(as.matrix(RNA_matrix[,intersect( basal_genes,rownames(brain_subset))]), na.rm = TRUE)
brain_subset
  gene_data_pos <- brain_subset@assays$SCT@data[intersect(classical_genes, rownames(brain_subset@assays$SCT@data)), ]
  gene_data_neg <- brain_subset@assays$SCT@data[intersect(basal_genes, rownames(brain_subset@assays$SCT@data)), ]
  
  brain_subset@meta.data$classical <- colMeans(gene_data_pos)
  brain_subset@meta.data$basal <- colMeans(gene_data_neg)

  img_name <- names(brain_subset@images)

  
  p_state <- SpatialDimPlot(
    brain_subset,
    images   = img_name,
    group.by = "state",
    cols     = cols_state,
    pt.size.factor = 1.6,
    label = FALSE
  ) + ggtitle(IMAGE)

  p_cluster <- SpatialDimPlot(
    brain_subset,
    images   = img_name,
    group.by = "seurat_clusters",
    #cols     = pal,
    pt.size.factor = 1.6,
    label = FALSE
  ) + ggtitle(IMAGE)
  
  p_hl <- SpatialDimPlot(
    brain_subset,
    images   = img_name,
    group.by = "is_cancer",
    cols     = c("TRUE"="#EEAEC2", "FALSE"="grey85"),
    pt.size.factor = 1.6,
    label = FALSE
  ) + ggtitle(paste0(IMAGE, " (Cancer highlight)"))
 
  SENE <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "sene",
    pt.size.factor = 1.6
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  
  EMT <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "emt",
    pt.size.factor = 1.6
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  #entropy <- SpatialFeaturePlot(
  #  brain_subset,
  #  images   = img_name,
  #  features = "entropy",
  #  pt.size.factor = 1.6
  #) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
  #  ggtitle(IMAGE)
  stemness <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "stemness",
    pt.size.factor = 1.6
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  basal <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "basal",
    pt.size.factor = 1.6
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  classical <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "classical",
    pt.size.factor = 1.6
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  
  immune_Supr <- SpatialFeaturePlot(
    brain_subset,
    images   = img_name,
    features = "immune_Supr",
    pt.size.factor = 1.6, max.cutoff=0.3
  ) +
    #scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA")) +
    ggtitle(IMAGE)
  
  plots_state[[IMAGE]]   <- p_state
  plots_cancerHL[[IMAGE]] <- p_hl
  plots_EMT[[IMAGE]] <- EMT
  plots_SENE[[IMAGE]] <- SENE
  plots_stemness[[IMAGE]] <- stemness
  #plots_entropy[[IMAGE]] <- entropy
  plots_basal[[IMAGE]] <- basal
  plots_classical[[IMAGE]] <- classical
  plots_cluster[[IMAGE]] <- p_cluster
  plots_immune[[IMAGE]] <- immune_Supr

  this_cells <- brain_subset@meta.data %>%
    as_tibble(rownames = "spot") %>%
    transmute(sample = IMAGE, spot, seurat_clusters,
              is_cancer = is_cancer, state = as.character(state),
              sene = sene, emt = emt, cancer = cancer, basal = basal, classical = classical,
              stemness = stemness)
  all_cells[[IMAGE]] <- this_cells
  
  lvl_state <- levels(brain_subset$state)
  this_comp <- this_cells %>%
    filter(is_cancer) %>%                            # only look at cancer cells
    mutate(state = factor(state, levels = lvl_state)) %>%
    count(sample, state, name = "cells") %>%
    group_by(sample) %>%
    complete(state = lvl_state, fill = list(cells = 0)) %>%
    mutate(prop = cells / sum(cells)) %>%
    ungroup() %>%
    filter(state != "Other")                          # no "Other" in cancer cells
  all_comp[[IMAGE]] <- this_comp
}

ggarrange(p_cluster,p_hl, SENE,EMT,immune_Supr,
          ncol = 5, nrow = 1,
          common.legend = TRUE, legend = "right"
)  

##1500*800 $$6, 9,13
image <- 'IU_PDA_T9'
p1 <-  plots_SENE[[image]]
p2 <-  plots_EMT[[image]]
p3 <- plots_stemness[[image]]
p4 <- plots_cancerHL[[image]]
p5 <- plots_classical[[image]]
p6 <- plots_basal[[image]]
p7 <- plots_cluster[[image]]
p7|p4|p1|p2|p3|p5|p6

df_cells <- bind_rows(all_cells)  # single-cell level details: one row per spot (this is the df you want)
df_comp  <- bind_rows(all_comp)   # proportion of four states per sample (slice) (cancer cells only)

library(patchwork )
N <- min(15, length(plots_state))
p_all <- wrap_plots(plots_stemness[c(1:15)], ncol = 5, guides = "collect") &
  theme(legend.position = "bottom", plot.title = element_text(size = 9))
print(p_all)


p_bar <- ggplot(df_comp, aes(x = sample, y = prop, fill = state)) +
  geom_col(width = 0.9, color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(values = cols_state, drop = FALSE, name = "State") +
  labs(x = "Sample", y = "Cell proportion (cancer cells)", title = "Cancer-cell states per sample") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
print(p_bar)

## 9) Self-check: whether multiple states exist within clusters (should be allowed; if all same, indicates you accidentally made "cluster-level labels")

# Only look at cancer cell boxplot (remove this line if you want to include all spots)
df_plot <- df_cells %>% filter(state != "Other",
                         !is.na(sample), !is.na(sene), !is.na(emt))

df_plot <-df_cells[df_cells$is_cancer ==TRUE, ]

# Long format
df_long <- df_plot %>%
  select(sample, sene, emt) %>%
  pivot_longer(c(sene, emt), names_to = "score", values_to = "value")

# Sample ordering: by overall median of both scores from high to low
ord <- df_long %>% group_by(sample) %>%
  summarise(med = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>% pull(sample)

ord <-c('IU_PDA_T2','IU_PDA_HM2','IU_PDA_HM2_2','IU_PDA_T3','IU_PDA_HM3',
        'IU_PDA_T4','IU_PDA_HM4','IU_PDA_T6','IU_PDA_LNM6','IU_PDA_HM6',
        'IU_PDA_T8','IU_PDA_LNM8','IU_PDA_HM8', 'IU_PDA_T9','IU_PDA_HM9',
        'IU_PDA_T10','IU_PDA_LNM10','IU_PDA_HM10','IU_PDA_T11','IU_PDA_HM11',
        'IU_PDA_T12','IU_PDA_HM12')
df_long <- df_long[df_long$sample %in% ord]
df_long$sample <- factor(df_long$sample, levels = ord)

# Plot
ggplot(df_long, aes(x = sample, y = value, fill = score)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3,
               position = position_dodge(width = 0.7)) +
  geom_jitter(aes(color = score),
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.7),
              size = 0.25, alpha = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = c(sene = "#F6A6C1", emt = "#86CFBA"),
                    name = NULL, labels = c("SENE", "EMT")) +
  scale_color_manual(values = c(sene = "#F6A6C1", emt = "#86CFBA")) +
  labs(x = "Sample", y = "Score",
       title = "Per-sample SENE & EMT scores (cancer spots)") +
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(
    axis.line   = element_line(color = "black", linewidth = 0.8),  
    axis.ticks  = element_line(linewidth = 1.0),                  
    axis.text   = element_text(size = 13, family = "Arial"),       
    axis.title  = element_text(size = 15, family = "Arial"),        
    plot.title  = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title= element_text(size = 13, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 8,family = "Arial"),
    axis.title.x = element_blank()
  )

df_plot$sample <- factor(df_plot$sample, levels = ord)
ggplot(df_plot, aes(x = sample, y = immune)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3) +
  geom_jitter(
    width = 0.15, height = 0,   # regular jitter
    size = 0.25, alpha = 0.25,
    color = "#F6A6C1", show.legend = FALSE
  )+
  labs(x = "Sample", y = "Score",
       title = "Per-sample SENE & EMT scores (cancer spots)") +
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(
    axis.line   = element_line(color = "black", linewidth = 0.8),  
    axis.ticks  = element_line(linewidth = 1.0),                  
    axis.text   = element_text(size = 13, family = "Arial"),       
    axis.title  = element_text(size = 15, family = "Arial"),        
    plot.title  = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title= element_text(size = 13, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 8,family = "Arial"),
    axis.title.x = element_blank()
  )


library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

## 1) Only use cancer cell spots; remove cluster_hl filter if you want all spots
df_use <- df_plot %>%
  filter(state != "Other",
         !is.na(sample), !is.na(sene), !is.na(emt)) %>%
  mutate(
    sample_short = sub("^IU_PDA_", "", sample),             # e.g. "T2", "HM2", "HM2_2"
    patient_id   = sub(".*?(\\d+).*", "\\1", sample_short), # extract number: "2"
    is_T_sample  = grepl("^T\\d+(?:$|_)", sample_short)     # whether it's a T sample
  )

## 2) Long format: SENE vs EMT
df_long1 <- df_use %>%
  select(sample, sample_short, patient_id, is_T_sample, sene, emt) %>%
  pivot_longer(c(sene, emt), names_to = "score", values_to = "value")

## 3) For each patient, each score: T_baseline vs other samples from that patient, pairwise Wilcoxon test
library(purrr)
library(stringr)
p_tab <- df_long1 %>%
  group_by(patient_id, score) %>%
  group_modify(function(dat, key){
    
    # Baseline T sample (if multiple T samples, only take the first; you can modify to select specific T)
    t_samps <- unique(dat$sample[dat$is_T_sample])
    if(length(t_samps) == 0) return(tibble())          # no T, skip
    t_base  <- t_samps[1]
    
    others  <- setdiff(unique(dat$sample), t_base)
    if(length(others) == 0) return(tibble())
    
    map_dfr(others, function(s2){
      x <- dat$value[dat$sample == t_base]
      y <- dat$value[dat$sample == s2]
      if(length(x) > 0 && length(y) > 0){
        wt <- suppressWarnings(wilcox.test(x, y))
        tibble(
          baseline   = t_base,
          compare    = s2,
          n_base     = length(x),
          n_comp     = length(y),
          median_base= median(x, na.rm = TRUE),
          median_comp= median(y, na.rm = TRUE),
          p          = wt$p.value
        )
      } else tibble()
    })
  }) %>%
  ungroup() %>%
  group_by(patient_id, score) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(score, patient_id, baseline, compare)

p_tab


library(patchwork)

# 1) Your defined display order (names must match names(plots_state) exactly)
img_order <- ord  # example, modify as needed

# 2) Only keep currently existing plots, maintain your defined order
keep <- img_order[img_order %in% names(plots_state)]
stopifnot(length(keep) > 0)

# (Optional) Only plot first N images
N <- 10
keep <- keep[11:20]

# (Optional) Add title=name to each plot
plots_ordered <- Map(`+`, plots_state[keep], lapply(keep, ggtitle))

# 3) Combine and theme
p_all <- (wrap_plots(plots_ordered, ncol = 11, guides = "collect") &
            theme(legend.position = "bottom",
                  plot.title = element_text(size = 9)))

# Display
p_all  # interactive environment direct display; or print(p_all)




library(Seurat)
library(dplyr)

genes <- c("EPCAM", "VIM")

m <- GetAssayData(brain_subset, slot = "data")[genes, ]  # log1p normalized

# Z-score genes by cells (rows: genes; columns: cells)
scale_mat <- t(scale(t(as.matrix(m))))  # standardize each gene row vector
E <- as.numeric(scale_mat["EPCAM", ])
M <- as.numeric(scale_mat["VIM", ])

# EMT direction
EMT_score <- M - E

# Plasticity (recommended: harmonic mean + small epsilon to prevent 0)
eps <- 1e-6
Plasticity_HM <- 2 / (1/(pmax(E,0)+eps) + 1/(pmax(M,0)+eps))

# Or use quantile to 0-1 then BPP
to01 <- function(x) (rank(x, ties.method="average")-1)/(length(x)-1)
Eq <- to01(E); Mq <- to01(M)
Plasticity_BPP <- 4 * Eq * Mq

#obj$EMT_score <- EMT_score
brain_subset$Plasticity <- as.numeric(Plasticity_BPP)  # or Plasticity_HM

# Visualization (UMAP)
FeaturePlot(brain_subset, features = c("Plasticity"), cols = c("navy","white","firebrick"))
stemness <- SpatialFeaturePlot(
  brain_subset,
  images   = img_name,
  features = "Plasticity",
  pt.size.factor = 1.6,
  max.cutoff=2

)+
  scale_fill_gradientn(colors = c("#F0F0F0", "#d0F0F0","#16CFBA"))


#####Mechanism 
####################below ###################
# Major cell type annotations for tumor/immune/stroma
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
TLS_genes  <- c("CXCL13","CCL19","CCL21","LTB","MS4A1","CD79A","ICOS","CXCR5")   # TLS representative gene set
cDC1_genes <- markers_list$cDC1
library(Seurat)
library(patchwork)

# 1) Select IMAGE name to plot (if you know it, directly write IMAGE <- "your_image_name")
IMAGE <- names(brain_subset@images)[1]

vln_by_cluster <- function(obj, features, group.by = "seurat_clusters", assay = "SCT",
                           ncol = 4, same.y.lims = TRUE, pt.size = 0.0) {
  feats <- intersect(features, rownames(obj))
  if (length(feats) == 0) stop("These genes are not in the object: ", paste(features, collapse=", "))
  if (length(feats) < length(features)) {
    message("The following genes were not found and will be skipped: ",
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

SpatialDimPlot(brain_subset,group.by="seurat_clusters")

cancer_clusters <- c("4","2","3") #HM2
cancer_clusters <- c("2") #HM2_2
cancer_clusters <- c("8","3") #T2
### cluster. ########
all_clusters <- levels(brain_subset$seurat_clusters)
cluster_colors <- setNames(rep("lightgray", length(all_clusters)), all_clusters)
cluster_colors[cancer_clusters] <- "#F5A5BE"
#cluster_colors[10] <- "blue"

c <- SpatialDimPlot(
  brain_subset,
  group.by = "seurat_clusters",
  images = IMAGE,
  cols = cluster_colors
)

# Panel genes (can be replaced with your library)
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

add_score <- function(obj, gs, nm) AddModuleScore(obj, features=list(intersect(gs, rownames(obj))), name=paste0(nm,"_"))
brain_subset <- add_score(brain_subset, CAF_barrier, "CAF")
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

# Spatial plots: primary vs metastasis, one page each (example)

d <- SpatialFeaturePlot(brain_subset, features=c("sene","emt","CAF_1","CXCL12ax_1","SASP_1","TGFb_1","cDC1_1","TLS_1"), pt.size.factor=1.6)


####################above ###################
SpatialFeaturePlot(brain_subset, features=c("_1"), pt.size.factor=1.6)
library(Seurat)
DefaultAssay(brain_subset) <- "SCT"

for (nm in names(markers_list)) {
  brain_subset <- AddModuleScore(brain_subset, features = list(markers_list[[nm]]), name = paste0(nm,"_MS_"))
}












ms_cols <- paste0(names(markers_list), "_MS_1")
ms_mat  <- FetchData(brain_subset, vars = ms_cols)
SpatialDimPlot(brain_subset, group.by = "celltype_anno", label = TRUE, pt.size.factor = 1.6)

call <- apply(ms_mat, 1, function(x){
  labs <- names(markers_list)
  labs[ which.max(x) ]
})
brain_subset$celltype_anno <- factor(call, levels = names(markers_list))

# Optional: binary classification of "Tumor" (tumor-enriched regions) for subsequent statistics
q75 <- quantile(brain_subset$Tumor_MS_1, 0.75, na.rm = TRUE)
brain_subset$tumor_region <- ifelse(brain_subset$Tumor_MS_1 >= q75, "TumorHigh","Others")




