
barplot(1,2)

rm(list=ls())
gc()
set.seed(1111)
setwd("/fs/ess/PCON0022/Weidong/Short/XY")


library(CellChat)
library(Matrix)
library(data.table) 
library(Cairo)


###################### P,p3 ###############################
cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_p_p3.rds")
  #TAM->TUMOR

#target_cell_types <- c("Zero","1","2","3","6","7","9","10","11","13","15","18")
  

source_list = c('1','14','9','17','3','5')
target_list = c('2','4','8','12','13','16')
interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(
  "TGFB1 - (TGFBR1+TGFBR2)",
  "TGFB1 - (ACVR1+TGFBR1)",
  "MIF - (CD74+CD44)",
  "IGF1 - (CD74+CD44)",
  "AREG - (EGFR+ERBB2)",
  "AREG - EGFR",
  "TNF - TNFRSF1A",
  "GZMA - PARD3",
  "GZMA - F2RL1",
  "IFNG - (IFNGR1+IFNGR2)"
))

# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.35, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p


gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks   # ggplot2 ≥ 3.3.0
# Filter the interaction_order to keep only those combinations that exist in x_levels
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
p = p + scale_x_discrete(limits = interaction_order_filtered)

#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(6.5,12.5
  #(i14_4 + i14_8)/2,
  #(i17_4 + i9_8)/2#,
  #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
               linewidth = 0.5, linetype = "dashed", color = "black") 



ggplot2::ggsave("p_p3.pdf", p, width = 5, height = 2, dpi = 300)



################ p p1 ################ 
cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_p_p1.rds")
#TAM->TUMOR
TAM<- c("5")
#target_cell_types <- c("Zero","1","2","3","6","7","9","10","11","13","15","18")

NK<- c("14") 
T8cell<- c("4") 
TAM<- c("5","8") 
B<- c("19","20") 

source_list = c('14','4','5','8','19','20')
target_list <- c("9","6","Zero","11","13","2","3","18","1","7","10","15")

interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(

  "MIF - (CD74+CD44)",
  "IGF1 - (CD74+CD44)",
  "AREG - (EGFR+ERBB2)",
  "AREG - EGFR",
  "TNF - TNFRSF1A",
  "GZMA - PARD3",
  "GZMA - F2RL1",
  "IFNG - (IFNGR1+IFNGR2)"
))

# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,#thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.35, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p
gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks
y_levels <- gb$layout$panel_params[[1]]$y$breaks

# Keep only the ones present in the plot
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
keep_rows_filtered <- keep_rows[keep_rows %in% y_levels]

remove_col = c('14 -> 13',
               '14 -> 6',
               '14 -> Zero',
               '19 -> 18',
               '20 -> 18',
               '19 -> 18',
               '4 -> 18',
               '5 -> 11',
               '5 -> 13',
               '5 -> 15',
               '5 -> 18',
               '5 -> 2',
               '5 -> 3',
               '5 -> 6',
               '5 -> 7',
               '8 -> 11',
               '8 -> 18',
               '8 -> 6',
               '8 -> Zero',
               '14 -> 11',
               '14 -> 18',
               '5 -> Zero',
               '8 -> 3'
               )
interaction_order_filtered <- setdiff(interaction_order_filtered, remove_col)

# Re-apply filters to both axes
p <- p +
  scale_x_discrete(limits = interaction_order_filtered, drop = TRUE) 

p
#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(6.5,12.5,20.5
            #(i14_4 + i14_8)/2,
            #(i17_4 + i9_8)/2#,
            #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
                   linewidth = 0.5, linetype = "dashed", color = "black") 
p


ggplot2::ggsave("p_p1.pdf", p, width = 6.5, height = 2, dpi = 300)


################ m p1  ##################

cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_m_p1.rds")
NK<- c("14") 
T8cell<- c("4") 
TAM<- c("5","8") 
B<- c("19","20") 
source_list = c('14','4','5','8','19','20')
target_list <- c("9","6","0","11","13","2","3","18","1","7","10","15")


interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(
  "MIF - (CD74+CD44)",
  "AREG - (EGFR+ERBB2)",
  "AREG - EGFR",
  "GZMA - PARD3",
  "GZMA - F2RL1"
))

# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.30, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p
gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks
y_levels <- gb$layout$panel_params[[1]]$y$breaks

# Keep only the ones present in the plot
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
keep_rows_filtered <- keep_rows[keep_rows %in% y_levels]

remove_col = c('14 -> 11',
               '14 -> 13',
               '14 -> 18',
               '19 -> 1',
               '19 -> 10',
               '19 -> 15',
               '19 -> 3',
               '19 -> 7',
               '20 -> 15',
               '20 -> 18',
               '4 -> 18',
               '5 -> 11',
               '5 -> 6',
               '5 -> 15',
               '5 -> 18',
               '8 -> 13',
               '8 -> 15',
               '8 -> 11',
               '8 -> 18',
               '14 -> 6',
               '8 -> 6'
              
)
interaction_order_filtered <- setdiff(interaction_order_filtered, remove_col)

# Re-apply filters to both axes
p <- p +
  scale_x_discrete(limits = interaction_order_filtered, drop = TRUE) 

p
#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(6.5,12.5
            #(i14_4 + i14_8)/2,
            #(i17_4 + i9_8)/2#,
            #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
                   linewidth = 0.5, linetype = "dashed", color = "black") 
p

ggplot2::ggsave("m_p1.pdf", p, width = 5, height = 1.8, dpi = 300)


################## p p2 ######################
cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_p_p2.rds")


NK<- c("4","20") 
T4cell<- c("14")
T8cell<- c("2","5","8","10")
TAM<- c("24") 
B<- c("18","21") 

target_cell_types <- c("11","22")

source_list = c('4','20','14','2','5','8','10','24','18','21')
target_list <- c("11","22")

interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(
  "MIF - (CD74+CD44)",
  "AREG - (EGFR+ERBB2)",
  "AREG - EGFR",
  "GZMA - PARD3",
  "GZMA - F2RL1"
))



# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.30, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p
gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks
y_levels <- gb$layout$panel_params[[1]]$y$breaks

# Keep only the ones present in the plot
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
keep_rows_filtered <- keep_rows[keep_rows %in% y_levels]

remove_col = c('14 -> 11'
               
)
interaction_order_filtered <- setdiff(interaction_order_filtered, remove_col)

# Re-apply filters to both axes
p <- p +
  scale_x_discrete(limits = interaction_order_filtered, drop = TRUE) 

p
#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(2.5,6.5
            #(i14_4 + i14_8)/2,
            #(i17_4 + i9_8)/2#,
            #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
                   linewidth = 0.5, linetype = "dashed", color = "black") 
p

ggplot2::ggsave("p_p2.pdf", p, width = 4, height = 1.8, dpi = 300)



################## m p2 #####################
NK<- c("4","20") 
T4cell<- c("14")
T8cell<- c("2","5","8","10")
TAM<- c("24") 
B<- c("18","21") 
target_cell_types <- c("11","22")

cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_m_p2.rds")



source_list = c('4','20','14','2','5','8','10','24','18','21')
target_list <- c("11","22")

interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(
  "MIF - (CD74+CD44)",
  "AREG - EGFR",
  "TNF - TNFRSF1A",
  "IFNG - (IFNGR1+IFNGR2)"
))


# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.30, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p
gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks
y_levels <- gb$layout$panel_params[[1]]$y$breaks

# Keep only the ones present in the plot
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
keep_rows_filtered <- keep_rows[keep_rows %in% y_levels]

remove_col = c('10 -> 11','10 -> 22','14 -> 11',
               '14 -> 22','18 -> 11','2 -> 11','20 -> 11','20 -> 22',
               '21 -> 11','5 -> 11','8 -> 11','5 -> 22','8 -> 22'
               
)
interaction_order_filtered <- setdiff(interaction_order_filtered, remove_col)

# Re-apply filters to both axes
p <- p +
  scale_x_discrete(limits = interaction_order_filtered, drop = TRUE) 

p
#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(2.5,5.5
            #(i14_4 + i14_8)/2,
            #(i17_4 + i9_8)/2#,
            #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
                   linewidth = 0.5, linetype = "dashed", color = "black") 
p

ggplot2::ggsave("m_p2.pdf", p, width = 4, height = 1.8, dpi = 300)



################## m p3 #####################

NK <- c("1","14")
T8cell <- c("3","5","7")
#T4cell<- c("9","17")
TAM<- c("0","19") 
#target_cell_types order : 8,13,16,2,12,4
source_list = c('1','14','9','17','3','5','7','Zero','19')
target_list = c('8','13','16','2','12','4')




cellchat<-readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchat_m_p3.rds")
#TAM->TUMOR

#target_cell_types <- c("Zero","1","2","3","6","7","9","10","11","13","15","18")

interaction_order <- unlist(lapply(source_list, function(source) {
  paste0(source, " -> ", target_list)
}))

# Reverse the order of keep_rows
keep_rows <- rev(c(
  "TGFB1 - (TGFBR1+TGFBR2)",
  "TGFB1 - (ACVR1+TGFBR1)",
  "MIF - (CD74+CD44)",
  
  "AREG - (EGFR+ERBB2)",
  "AREG - EGFR",
  "TNF - TNFRSF1A",
  "GZMA - PARD3",
  "GZMA - F2RL1"

))

# Apply this reversed order to the plot
p <- netVisual_bubble(cellchat, sources.use = source_list,remove.isolate = TRUE,
                      grid.on = FALSE,thresh = 0.05,
                      targets.use = target_list, return.data = FALSE) +
  scale_color_gradientn(colours = c("#FADB7D", "#FF3FB7"),
                        name = "Commun. Prob.") +
  theme(panel.grid = element_blank(),axis.text.y = element_text(size = 4),
        legend.key.height = unit(0.35, "cm"),
        axis.text.x = element_text(size = 4)) +
  scale_y_discrete(limits = keep_rows, drop = TRUE) 
p


gb <- ggplot_build(p)
x_levels <- gb$layout$panel_params[[1]]$x$breaks
y_levels <- gb$layout$panel_params[[1]]$y$breaks

# Keep only the ones present in the plot
interaction_order_filtered <- interaction_order[interaction_order %in% x_levels]
keep_rows_filtered <- keep_rows[keep_rows %in% y_levels]

remove_col = c('1 -> 13','1 -> 16','17 -> 12','1 -> 8','9 -> 8',
               '17 -> 4','17 -> 8','19 -> 12','19 -> 13','19 -> 16',
               '3 -> 13','3 -> 16','5 -> 13','5 -> 16','5 -> 8','7 -> 4','Zero -> 13','Zero -> 16','Zero -> 8'
               
)
interaction_order_filtered <- setdiff(interaction_order_filtered, remove_col)

# Re-apply filters to both axes
p <- p +
  scale_x_discrete(limits = interaction_order_filtered, drop = TRUE) 

#i14_4 <- which(x_levels == "3 -> 2")
#i14_8 <- which(x_levels == "17 -> 12")
#i9_8  <- which(x_levels == "14 -> 12")
#i17_4 <- which(x_levels == "9 -> 2")
#i3_8  <- which(x_levels == "3 -> 8")
#i5_12  <- which(x_levels == "5 -> 12")

# compute midpoints
midpts <- c(2.5,10.5,6.5
            #(i14_4 + i14_8)/2,
            #(i17_4 + i9_8)/2#,
            #(i5_12 + i3_8)/2
)

p = p + geom_vline(xintercept = midpts,
                   linewidth = 0.5, linetype = "dashed", color = "black") 
p


ggplot2::ggsave("m_p3.pdf", p, width = 5, height = 2, dpi = 300)




######################### p1 #################

library(RColorBrewer)
# Choose as many colors as tumor clusters

setwd('/fs/ess/PCON0022/Weidong/Short/XY/')
H_tumor <- c('Tumor1','Tumor2','Tumor3','Tumor7','Tumor10','Tumor13','Tumor15')
L_tumor <- c('Tumor6','Tumor0','Tumor9','Tumor11')
cellchat = readRDS('/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp1_full.rds')

tumor_clusters <- grep("^Tumor", levels(cellchat@idents), value = TRUE)
#cluster.colors <- setNames(colorRampPalette(my.colors)(length(tumor_clusters)), pretty_names)

n_tumor <- length(tumor_clusters)
my.colors <- brewer.pal(n = max(3, min(12, n_tumor)), name = "Pastel1")

# Make a named vector so each cluster always gets the same color

cluster.colors <- setNames(colorRampPalette(my.colors)(n_tumor), tumor_clusters)

gene_list = c("EGF","CypA","PARs","SPP1","APP","LAMININ","GALECTIN","MK")
pdf('/fs/ess/PCON0022/Weidong/Short/XY/CCC_figure/aggregate/p1.pdf', 
    width = 6, height = 32)   # adjust width/height for readability

# set layout: 1 column, 8 rows
par(mfrow = c(8, 1), 
    mar = c(1, 1, 2, 1),   # bottom, left, top, right
    oma = c(0, 0, 0, 0),   # outer margins
    mgp = c(2, 0.5, 0))    # axis title/label/line spacing

for (gene in gene_list[1:8]) {  # first 8 genes, change if needed
  netVisual_aggregate(cellchat, 
                      signaling = gene, 
                      layout = "chord", 
                      sources.use = tumor_clusters, 
                      targets.use = tumor_clusters,
                      remove.isolate = TRUE,
                      color.use = cluster.colors,
                      vertex.label.cex = 2.5,
                      vertex.label = pretty_names)
  #title(main = paste("Signaling:", gene), line = -1)  # optional label
}

dev.off()


#####################
library(RColorBrewer)

cellchat1 <- readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp1_full.rds")
cellchat2 <- readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp2_full.rds")
cellchat <- readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp3_full.rds")
####p2
e <- netVisual_bubble(cellchat2, sources.use = "Tumor11", 
                      targets.use = "Tumor12",return.data = TRUE)
f <- netVisual_bubble(cellchat2, sources.use = "Tumor12", 
                      targets.use = "Tumor11",return.data = TRUE)

####p1



####p3

H_tumor <- c('Tumor2','Tumor4','Tumor12')
L_tumor <- c('Tumor13','Tumor8','Tumor16')

a <- netVisual_bubble(cellchat, sources.use = H_tumor, 
                      targets.use = L_tumor,return.data = TRUE)

b<-netVisual_bubble(cellchat, sources.use = L_tumor, 
                    targets.use = H_tumor,return.data = TRUE)


unique(setdiff(a$communication$pathway_name,b$communication$pathway_name))
unique(setdiff(b$communication$pathway_name,a$communication$pathway_name))




cellchat <- readRDS("/fs/ess/PAS1475/Xiaoying/Pancreas_meta/P1/cellchatp1_full.rds")
H_tumor <- c('Tumor1','Tumor2','Tumor3','Tumor7','Tumor10','Tumor13','Tumor15')
L_tumor <- c('Tumor6','Tumor0','Tumor9','Tumor11')

c<-netVisual_bubble(cellchat1, sources.use = H_tumor, 
                    targets.use = L_tumor,return.data = TRUE)
d <-netVisual_bubble(cellchat1, sources.use = L_tumor, 
                     targets.use = H_tumor,return.data = TRUE)

id_p = 'p1cd'
all_list = unique(setdiff(c$communication$pathway_name,d$communication$pathway_name))
tumor_clusters <- grep("^Tumor", levels(cellchat@idents), value = TRUE)
#cluster.colors <- setNames(colorRampPalette(my.colors)(length(tumor_clusters)), pretty_names)

n_tumor <- length(tumor_clusters)
my.colors <- brewer.pal(n = max(3, min(12, n_tumor)), name = "Pastel1")

# Make a named vector so each cluster always gets the same color

cluster.colors <- setNames(colorRampPalette(my.colors)(n_tumor), tumor_clusters)

gene_list = c("EGF","CypA","PARs","SPP1","APP","LAMININ","GALECTIN","MK")
gene_list <- na.omit(setdiff(all_list, gene_list))
pdf(paste0('/fs/ess/PCON0022/Weidong/Short/XY/CCC_figure/sup_aggregate/',id_p,'.pdf'), 
    width = 6, height = 32)   # adjust width/height for readability

# set layout: 1 column, 8 rows
par(mfrow = c(8, 1), 
    mar = c(1, 1, 2, 1),   # bottom, left, top, right
    oma = c(0, 0, 0, 0),   # outer margins
    mgp = c(2, 0.5, 0))    # axis title/label/line spacing

for (gene in gene_list) {  # first 8 genes, change if needed
  netVisual_aggregate(cellchat, 
                      signaling = gene, 
                      layout = "chord", 
                      sources.use = tumor_clusters, 
                      targets.use = tumor_clusters,
                      remove.isolate = TRUE,
                      color.use = cluster.colors,
                      vertex.label.cex = 1.5)
  #title(main = paste("Signaling:", gene), line = -1)  # optional label
}

dev.off()





###################
