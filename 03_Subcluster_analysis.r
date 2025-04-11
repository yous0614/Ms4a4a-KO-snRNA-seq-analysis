conda activate seurat_v5
R

library(qs)
library(Seurat)
library(dplyr)
library(ggrepel)
library(stringr)
library(Azimuth)

#########################################################
# load the data
#########################################################

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types")
Astro = qread("MS4A_KO_mice_Astrocytes_HarmonyIntegrated.qs"        )
Excit = qread("MS4A_KO_mice_Excitatory_neurons_HarmonyIntegrated.qs")
Inhib = qread("MS4A_KO_mice_Inhibitory_neurons_HarmonyIntegrated.qs")
Micro = qread("MS4A_KO_mice_Microglia_HarmonyIntegrated.qs"         )
Oligo = qread("MS4A_KO_mice_Oligodendrocyte_HarmonyIntegrated.qs"   )
OPC = qread("MS4A_KO_mice_OPC_HarmonyIntegrated.qs"               )
Other = qread("MS4A_KO_mice_Other_HarmonyIntegrated.qs" )

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis")
cell_types = c("0-Microglia", "1-Astrocytes", "2-Excitatory_neurons", "3-Inhibitory_neurons", "4-Oligodendrocyte", "5-OPC", "6-Other")
for(cell in cell_types){
  dir.create(file.path(cell, "0-UMAP"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "1-Proportion"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "2-FindAllMarkers"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "3-DEG_analysis"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "4-Volcano_plot"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "5-EnrichR"), showWarnings = FALSE, recursive = TRUE)
}

#########################################################
# UMAP
#########################################################

obj_list = list(Micro, Astro, Excit, Inhib, Oligo, OPC, Other)
for(i in 1:length(obj_list)){
  obj = obj_list[[i]]
  cell = cell_types[i]
  png(paste0(cell, "/0-UMAP/", cell, "_UMAP.png"), width = 6, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7) + 
    labs(title = paste0(cell, " Subclusters")) +
    theme(text = element_text(size=20, face="bold"))
  print(p)
  dev.off()
}

#########################################################
# Proportion
#########################################################

library(dplyr)

# Assume obj_list is a named list: list("Microglia" = micro_obj, "Astrocytes" = astro_obj, ...)

for(i in 1:length(obj_list)){

  obj = obj_list[[i]]
  cell = cell_types[i]

  output_dir <- paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/1-Proportion")
  df <- obj@meta.data
  
  # Count cell types per animal
  df_summary <- df %>%
    group_by(orig.ident, harmony_clusters, genotype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(percentage = n / sum(n)) %>%
    ungroup()
  
  # Save CSV
  write.csv(df_summary, file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_proportion.csv")), row.names = FALSE)
  
  # Do pairwise t-tests for each cell type
  cat("\n==========", cell, "==========\n")
  
  for (j in unique(df_summary$harmony_clusters)){
    test <- df_summary %>% filter(harmony_clusters == j)
    if (length(unique(test$genotype)) > 1) {  # Ensure both genotypes are present
      res <- pairwise.t.test(test$percentage, test$genotype, paired = FALSE, p.adjust.method = "none", pool.sd = FALSE)
      pval <- res$p.value
      cat("Cell type:", j, "| p-value:\n")
      print(pval)
    } else {
      cat("Cell type:", j, "| Only one genotype present — skipping test.\n")
    }
  }
}


# No significant cluster

#########################################################
# DEG on WT vs KO
#########################################################

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
protein_coding_genes <- getBM(attributes = c('mgi_symbol'), 
                              filters = 'biotype', 
                              values = 'protein_coding', 
                              mart = ensembl)
protein_coding_genes <- protein_coding_genes$mgi_symbol
head(protein_coding_genes)
 
Idents(Micro) <- "genotype"

Mic.0 = Micro %>% subset(harmony_clusters == 0) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.1 = Micro %>% subset(harmony_clusters == 1) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.2 = Micro %>% subset(harmony_clusters == 2) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.3 = Micro %>% subset(harmony_clusters == 3) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.4 = Micro %>% subset(harmony_clusters == 4) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.5 = Micro %>% subset(harmony_clusters == 5) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.6 = Micro %>% subset(harmony_clusters == 6) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.7 = Micro %>% subset(harmony_clusters == 7) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 
Mic.8 = Micro %>% subset(harmony_clusters == 8) %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) 

Mic.0 = Mic.0 %>% mutate(gene_biotype = ifelse(row.names(Mic.0) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.0))
Mic.1 = Mic.1 %>% mutate(gene_biotype = ifelse(row.names(Mic.1) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.1))
Mic.2 = Mic.2 %>% mutate(gene_biotype = ifelse(row.names(Mic.2) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.2))
Mic.3 = Mic.3 %>% mutate(gene_biotype = ifelse(row.names(Mic.3) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.3))
Mic.4 = Mic.4 %>% mutate(gene_biotype = ifelse(row.names(Mic.4) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.4))
Mic.5 = Mic.5 %>% mutate(gene_biotype = ifelse(row.names(Mic.5) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.5))
Mic.6 = Mic.6 %>% mutate(gene_biotype = ifelse(row.names(Mic.6) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.6))
Mic.7 = Mic.7 %>% mutate(gene_biotype = ifelse(row.names(Mic.7) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.7))
Mic.8 = Mic.8 %>% mutate(gene_biotype = ifelse(row.names(Mic.8) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.8)) 

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/0-Microglia/3-DEG_analysis")

Mic.0 %>% write.csv("2025-4-6 KO_vs_WT_Mic.0.csv")
Mic.1 %>% write.csv("2025-4-6 KO_vs_WT_Mic.1.csv")
Mic.2 %>% write.csv("2025-4-6 KO_vs_WT_Mic.2.csv")
Mic.3 %>% write.csv("2025-4-6 KO_vs_WT_Mic.3.csv")
Mic.4 %>% write.csv("2025-4-6 KO_vs_WT_Mic.4.csv")
Mic.5 %>% write.csv("2025-4-6 KO_vs_WT_Mic.5.csv")
Mic.6 %>% write.csv("2025-4-6 KO_vs_WT_Mic.6.csv")
Mic.7 %>% write.csv("2025-4-6 KO_vs_WT_Mic.7.csv")
Mic.8 %>% write.csv("2025-4-6 KO_vs_WT_Mic.8.csv")

Mic.0 = read.csv("2025-4-6 KO_vs_WT_Mic.0.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.1 = read.csv("2025-4-6 KO_vs_WT_Mic.1.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.2 = read.csv("2025-4-6 KO_vs_WT_Mic.2.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.3 = read.csv("2025-4-6 KO_vs_WT_Mic.3.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.4 = read.csv("2025-4-6 KO_vs_WT_Mic.4.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.5 = read.csv("2025-4-6 KO_vs_WT_Mic.5.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.6 = read.csv("2025-4-6 KO_vs_WT_Mic.6.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.7 = read.csv("2025-4-6 KO_vs_WT_Mic.7.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.8 = read.csv("2025-4-6 KO_vs_WT_Mic.8.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

# Function to process a single cell type
process_DEG_by_cluster <- function(seurat_obj, celltype_name, protein_coding_genes, output_dir) {
  Idents(seurat_obj) <- "genotype"
  
  # Get the number of unique harmony clusters
  clusters <- sort(unique(seurat_obj$harmony_clusters))
  
  for (clust in clusters) {
    # Perform differential expression
    markers <- seurat_obj %>%
      subset(harmony_clusters == clust) %>%
      FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE)
    
    # Add gene info and filter
    markers <- markers %>%
      mutate(
        gene = rownames(markers),
        gene_biotype = ifelse(rownames(markers) %in% protein_coding_genes, "protein_coding_genes", NA)
      )
    
    # Define filename
    out_csv <- file.path(output_dir, paste0("3-DEG_analysis/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    
    # Write to CSV
    write.csv(markers, out_csv)
    
    #Volcano plot
    output_vol <- file.path(output_dir, paste0("4-Volcano_plot/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    alternative_volcano_plot(markers, "p_val_adj", "avg_log2FC", path = output_vol, title = paste0(celltype_name, "_KO_vs_WT_", clust), genotype = "KO vs WT", celltype = celltype_name)

    #EnrichR
    output_path <- file.path(output_dir, paste0("5-EnrichR/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    enrichR_pathway(markers, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = output_path, title = paste0(celltype_name, "_KO_vs_WT_", clust))

    # Optional: Read back and filter (if continuing with further analysis)
    #assign(paste0(celltype_name, "_", clust), read.csv(out_csv, row.names = 1) %>% filter(gene_biotype == "protein_coding_genes"), envir = .GlobalEnv)
  }
}

process_DEG_by_cluster(Micro, "Microglia", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/0-Microglia")
process_DEG_by_cluster(Astro, "Astrocytes", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/1-Astrocytes")
process_DEG_by_cluster(Excit, "Excitatory_neurons", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/2-Excitatory_neurons")
process_DEG_by_cluster(Inhib, "Inhibitory_neurons", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/3-Inhibitory_neurons")
process_DEG_by_cluster(Oligo, "Oligodendrocytes", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/4-Oligodendrocyte")
process_DEG_by_cluster(OPC, "OPC", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/5-OPC")
process_DEG_by_cluster(Other, "Other", protein_coding_genes, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/6-Other")







source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/3-Volcano_plot/"

Color = hue_pal()(7)

volcano_plot(Ext.de, "p_val_adj", "avg_log2FC", path = Path, title = "Excitatory_neuron_KO_vs_WT", x = 25, y = 100)
alternative_volcano_plot(Ext.de, "p_val_adj", "avg_log2FC", path = Path, title = "Excitatory_neuron_KO_vs_WT", genotype = "KO vs WT", celltype = "Excitatory Neuron", xlimit = 10, color = Color[1])
alternative_volcano_plot(Inh.de, "p_val_adj", "avg_log2FC", path = Path, title = "Inhibitory_neuron_KO_vs_WT", genotype = "KO vs WT", celltype = "Inhibitory Neuron", color = Color[2])
alternative_volcano_plot(Ast.de, "p_val_adj", "avg_log2FC", path = Path, title = "Astrocyte_KO_vs_WT", genotype = "KO vs WT", celltype = "Astrocyte", color = Color[3])
alternative_volcano_plot(Mic.de, "p_val_adj", "avg_log2FC", path = Path, title = "Microglia_KO_vs_WT", genotype = "KO vs WT", celltype = "Microglia", color = Color[4])
alternative_volcano_plot(Oli.de, "p_val_adj", "avg_log2FC", path = Path, title = "Oligodendrocyte_KO_vs_WT", genotype = "KO vs WT", celltype = "Oligodendrocyte", color = Color[5])
alternative_volcano_plot(Opc.de, "p_val_adj", "avg_log2FC", path = Path, title = "OPC_KO_vs_WT", genotype = "KO vs WT", celltype = "OPC", color = Color[6])
alternative_volcano_plot(Oth.de, "p_val_adj", "avg_log2FC", path = Path, title = "Other_KO_vs_WT", genotype = "KO vs WT", celltype = "Other", color = Color[7])






#########################################################
# Identify microglia subclusters 
#########################################################

# Microglia
Micro = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Microglia_HarmonyIntegrated.qs")
Idents(Micro) <- "harmony_clusters"  # Set the identity class to harmony_clusters
FindAllMarkers(Micro, only.pos = TRUE) %>% filter(p_val_adj < 0.05) %>% write.csv("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/2-Subclusters/0-Microglia/0-UMAP/0-Microglia_FindAllMarkers.csv")

png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/2-Subclusters/0-Microglia/0-UMAP/Microglia_UMAP.png", width = 6, height = 5, units = "in", res = 300)
DimPlot(Micro, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7) + 
  labs(title = "Microglia Subclusters") +
  theme(text = element_text(size=20, face="bold"))
dev.off()
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/2-Subclusters/0-Microglia/0-UMAP/Microglia_UMAP_markers.png",  width = 25, height = 20, units = "in", res = 300)
FeaturePlot(Micro, features = c("Tmem119", "Aif1", "Cx3cr1", "P2ry12", "Trem2", "Cd68", "Clec7a", "Ms4a4a", "Ms4a7", "Cd163", "Apoe", "Ifit3"), reduction = "umap.harmony", order = TRUE, pt.size = 0.7, ncol = 4)
dev.off()
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/2-Subclusters/0-Microglia/0-UMAP/Microglia_UMAP_markers2.png",  width = 25, height = 20, units = "in", res = 300)
FeaturePlot(Micro, features = c("Tmem119", "Aif1", "Cx3cr1", "P2ry12", "Trem2", "Cd3e", "Cd3d", "Cd4", "Cd8a", "Ms4a4a", "Ms4a7", "Cd163"), reduction = "umap.harmony", order = TRUE, pt.size = 0.7, ncol = 4)
dev.off()

#Rename the clusters
# DF = Trem2PIC_microglia@meta.data
# DF <- DF %>%
#   mutate(microglia_subcluster = case_when(
#     harmony_clusters == 0 ~ "Mic.0",
#     harmony_clusters == 1 ~ "Mic.1",
#     harmony_clusters == 2 ~ "Mic.2",
#     harmony_clusters == 3 ~ "BAM",
#     harmony_clusters == 4 ~ "Mic.3",
#     TRUE ~ as.character(harmony_clusters)
#   ))
# Trem2PIC_microglia = AddMetaData(Trem2PIC_microglia, DF$microglia_subcluster, col.name = "microglia_subcluster")

#Marker for macrophage
Idents(Trem2PIC_microglia) = "microglia_subcluster"
Mac_markers = FindMarkers(Trem2PIC_microglia, ident.1 = "BAM", only.pos = T)
Mic_markers = FindAllMarkers(Trem2PIC_microglia %>% subset(microglia_subcluster != "BAM"), only.pos = T)
Mac_markers %>% write.csv("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2PIC_microglia_BAM_Markers.csv")
Mic_markers %>% write.csv("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2PIC_microglia_Microglia_Markers.csv")

# Set the desired order of identities
desired_order <- c("Mic.0", "Mic.1", "Mic.2", "Mic.3", "BAM")  # Example order
Trem2PIC_microglia@meta.data$microglia_subcluster <- factor(Trem2PIC_microglia@meta.data$microglia_subcluster, levels = desired_order)
levels(Idents(Trem2PIC_microglia))

png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_microglia_clusters.png", width = 6, height = 5, units = "in", res = 300)
DimPlot(Trem2PIC_microglia, group.by = "microglia_subcluster", pt.size = 0.7, reduction = "umap.harmony_20", label = F) + labs(title = "")
dev.off()

png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_manual_annotation.png", width = 7, height = 5, units = "in", res = 300)
DimPlot(obj, group.by = "manual_annotation", pt.size = 0.7, reduction = "umap.harmony_20", label = F) + labs(title = "")
dev.off()

Mac_markers = FindMarkers(Trem2PIC_microglia, ident.1 = "BAM", group.by = "microglia_subcluster", only.pos = T)
Mic_markers = FindAllMarkers(Trem2PIC_microglia %>% subset(microglia_subcluster != "BAM"), only.pos = T)
Mac_markers %>% write.csv("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2PIC_microglia_BAM_Markers.csv")
Mic_markers %>% write.csv("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2PIC_microglia_Microglia_Markers.csv")

#Generate heatmap for all the microglia subclusters

Trem2PIC_microglia <- ScaleData(Trem2PIC_microglia, features = rownames(Trem2PIC_microglia))
Mac_markers %>% filter(p_val_adj < 0.05) %>% head(20) %>% rownames() %>% unique -> BAM_top20
Mic_markers %>% filter(cluster == "Mic.0" & p_val_adj < 0.05)%>% head(20) %>% rownames() %>% unique -> Mic.0_top20
Mic_markers %>% filter(cluster == "Mic.1" & p_val_adj < 0.05)%>% head(20) %>% rownames() %>% unique -> Mic.1_top20
Mic_markers %>% filter(cluster == "Mic.2" & p_val_adj < 0.05)%>% head(20) %>% rownames() %>% unique -> Mic.2_top20
Mic_markers %>% filter(cluster == "Mic.3" & p_val_adj < 0.05)%>% head(20) %>% rownames() %>% unique -> Mic.3_top20

Heatmap_list = c(Mic.0_top20, Mic.1_top20, Mic.2_top20, Mic.3_top20, BAM_top20)
maxcells <- min(table(Idents(Trem2PIC_microglia)))

png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Heatmap for microglia subclusters_named_only.png" , width = 10, height = 9, units = "in", res = 300)
DoHeatmap(subset(Trem2PIC_microglia, downsample = maxcells), group.by = "microglia_subcluster", features = Heatmap_list)
dev.off()

qsave(Trem2PIC_microglia, "/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/0-data/0-Whole_brain/2025-3-1 microglia_BAM_integrated_20PCs_res0.2_SF.qs")

#########################################################
# Figure representation for each of the microglia subclusters
#########################################################

Trem2PIC_microglia = qread("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/0-data/0-Whole_brain/2025-3-1 microglia_BAM_integrated_20PCs_res0.2_SF.qs")

#Mic0
Mic0_marker = c("P2ry12", "Tmem119", "P2ry13", "Cx3cr1")
png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_Mic.0_homeostatic.png", width = 5, height = 5, units = "in", res = 300)
p = Trem2PIC_microglia %>% FeaturePlot(Mic0_marker, pt.size = 0.7, reduction = "umap.harmony_20", order = T, min.cutoff = c(3.5, 2.5, 2.5, 3.5), ncol = 2) & theme(text = element_text(size=20,face="bold"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
print(p)
dev.off()

#Mic1
Mic1_marker = c("Ifit3", "Ifitm3", "Stat1", "Irf7")
png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_Mic.1_interferon.png", width = 5, height = 5, units = "in", res = 300)
p = Trem2PIC_microglia %>% FeaturePlot(Mic1_marker, pt.size = 0.7, reduction = "umap.harmony_20", order = T, min.cutoff = c(0,1,2.1,0), ncol = 2) & theme(text = element_text(size=20,face="bold"))
for(i in 1:length(p)) {
p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
print(p)
dev.off()

#Mic2
Mic2_marker = c("Ccl3", "Cd63", "Ccl4", "Cd9")
png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_Mic.2_inflammatory.png", width = 5, height = 5, units = "in", res = 300)
p = Trem2PIC_microglia %>% FeaturePlot(Mic2_marker, pt.size = 0.7, reduction = "umap.harmony_20", order = T, min.cutoff = c(0,3,0,3), ncol = 2) & theme(text = element_text(size=20,face="bold"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
print(p)  # Increase the font size to 14
dev.off()


#Mic3
Mic3_marker = c("Robo1", "Robo2", "Gabrg2", "Gabrg3")
png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_Mic.3_axon_guidance.png", width = 5, height = 5, units = "in", res = 300)
p = Trem2PIC_microglia %>% FeaturePlot(Mic3_marker, pt.size = 0.7, reduction = "umap.harmony_20", order = T, min.cutoff = c(0,3,0,3), ncol = 2) & theme(text = element_text(size=20,face="bold"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
print(p)  # Increase the font size to 14
dev.off()


#BAM
BAM_marker = c("Cd163", "Mrc1", "Clec10a", "Ms4a7")
png("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/2025-3-1 Trem2_PolyIC_UMAP_ScaleData_harmony_20_res0.2_BAM.png", width = 5, height = 5, units = "in", res = 300)
p = Trem2PIC_microglia %>% FeaturePlot(BAM_marker, pt.size = 0.7, reduction = "umap.harmony_20", order = T, min.cutoff = c(0,3,0,0), ncol = 2) & theme(text = element_text(size=20,face="bold"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
print(p)  # Increase the font size to 14
dev.off()


