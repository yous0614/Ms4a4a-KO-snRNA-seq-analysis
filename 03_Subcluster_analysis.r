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
  dir.create(file.path(cell, "0-UMAP", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "0-UMAP", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "0-UMAP", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(cell, "1-Proportion", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "1-Proportion", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "1-Proportion", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)

  dir.create(file.path(cell, "2-FindAllMarkers", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "2-FindAllMarkers", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "2-FindAllMarkers", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(cell, "3-DEG_analysis", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "3-DEG_analysis", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "3-DEG_analysis", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(cell, "4-Volcano_plot", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "4-Volcano_plot", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "4-Volcano_plot", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(cell, "5-EnrichR", "All_cell_type", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "5-EnrichR", "All_cell_type", "Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(cell, "5-EnrichR", "Female_only", "5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT"), showWarnings = FALSE, recursive = TRUE)
}

obj_list = list(Micro, Astro, Excit, Inhib, Oligo, OPC, Other)

for(i in 1:length(obj_list)){
  merged_obj <- obj_list[[i]]
  merged_obj$fAD = ifelse(str_detect(merged_obj$orig.ident, "Exp12"), "Control", "5xFAD")
  merged_obj$condition = paste0(merged_obj$fAD, "_", merged_obj$genotype)
  merged_obj@meta.data %>% head
  merged_obj %>% qsave(file.path("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types", paste0("MS4A_KO_mice_Harmony_Integration_", cell_types[i], ".qs")))
}

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types")
Micro = qread("MS4A_KO_mice_Harmony_Integration_0-Microglia.qs"         )
Astro = qread("MS4A_KO_mice_Harmony_Integration_1-Astrocytes.qs"        )
Excit = qread("MS4A_KO_mice_Harmony_Integration_2-Excitatory_neurons.qs")
Inhib = qread("MS4A_KO_mice_Harmony_Integration_3-Inhibitory_neurons.qs")
Oligo = qread("MS4A_KO_mice_Harmony_Integration_4-Oligodendrocyte.qs"   )
OPC   = qread("MS4A_KO_mice_Harmony_Integration_5-OPC.qs"               )
Other = qread("MS4A_KO_mice_Harmony_Integration_6-Other.qs"             )

##################################################################################################################
##################################################################################################################
###                                                                                                           ####
###                                       All the genotype and sex                                            ####
###                                                                                                           ####
##################################################################################################################
##################################################################################################################

#---------------------------------------------------------------- 5xFAD_Ms4a4a KO vs 5xFAD_Ms4a4a WT ----------------------------------------------------------------

#########################################################
# UMAP
#########################################################

obj_list = list(Micro, Astro, Excit, Inhib, Oligo, OPC, Other)
setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/")

for(i in 1:length(obj_list)){
  obj = obj_list[[i]] %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
  cell = cell_types[i]
  png(paste0(cell, "/0-UMAP/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", cell, "_UMAP.png"), width = 6, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7) + 
    labs(title = paste0(cell, " Subclusters")) +
    theme(text = element_text(size=20, face="bold"))
  print(p)
  dev.off()
  png(paste0(cell, "/0-UMAP/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", cell, "_UMAP_KO_vs_WT.png"), width = 12, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7, split.by = "condition") + 
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

  obj = obj_list[[i]] %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
  cell = cell_types[i]

  output_dir <- paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/1-Proportion/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/")
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
# FindAllMarkers
#########################################################

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

for(i in 1:length(obj_list)){
  obj = obj_list[[i]] %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
  cell = cell_types[i]

  Idents(obj) <- "harmony_clusters"
  DEG_df = FindAllMarkers(obj, only.pos = TRUE) 
  DEG_df %>% write.csv(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", Sys.Date(), "_", cell, "_FindAllMarkers.csv"))
  
  for(j in 1:length(unique(DEG_df$cluster))){
    clust = unique(DEG_df$cluster)[j]
    output_path <- file.path(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"))
    enrichR_pathway(DEG_df %>% filter(cluster == clust), pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = output_path, title = paste0(cell, "_", clust))
  }
}


#########################################################
# DEG analysis on 5xFAD_Ms4a4a_KO vs 5xFAD_Ms4a4a_WT
#########################################################

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
protein_coding_genes <- getBM(attributes = c('mgi_symbol'), 
                              filters = 'biotype', 
                              values = 'protein_coding', 
                              mart = ensembl)
protein_coding_genes <- protein_coding_genes$mgi_symbol
head(protein_coding_genes)
 
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

# Function to process a single cell type
process_DEG_by_cluster <- function(seurat_obj, celltype_name, protein_coding_genes, output_dir) {
  seurat_obj = seurat_obj %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
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
    out_csv <- file.path(output_dir, paste0("3-DEG_analysis/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    
    # Write to CSV
    write.csv(markers, out_csv)
    
    #Volcano plot
    output_vol <- file.path(output_dir, paste0("4-Volcano_plot/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    alternative_volcano_plot(markers, "p_val_adj", "avg_log2FC", path = output_vol, title = paste0(celltype_name, "_KO_vs_WT_", clust), genotype = "KO vs WT", celltype = celltype_name)

    #EnrichR
    output_path <- file.path(output_dir, paste0("5-EnrichR/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
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

#---------------------------------------------------------------- Control_Ms4a4a KO vs Control_Ms4a4a WT ----------------------------------------------------------------

#########################################################
# UMAP
#########################################################

obj_list = list(Micro, Astro, Excit, Inhib, Oligo, OPC, Other)
setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/")

for(i in 1:length(obj_list)){
  obj = obj_list[[i]] %>% subset(condition == "Control_KO" | condition == "Control_WT")
  cell = cell_types[i]
  png(paste0(cell, "/0-UMAP/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", cell, "_UMAP.png"), width = 6, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7) + 
    labs(title = paste0(cell, " Subclusters")) +
    theme(text = element_text(size=20, face="bold"))
  print(p)
  dev.off()
  png(paste0(cell, "/0-UMAP/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", cell, "_UMAP_KO_vs_WT.png"), width = 12, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7, split.by = "condition") + 
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
  tryCatch({
    obj = obj_list[[i]] %>% subset(condition == "Control_KO" | condition == "Control_WT")
    cell = cell_types[i]

    output_dir <- paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/1-Proportion/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/")
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
  }, error = function(e) {
    cat("An error occurred:", conditionMessage(e), "\n")
  })
}


# No significant cluster


#########################################################
# FindAllMarkers
#########################################################

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

for(i in 1:length(obj_list)){
  tryCatch({
  obj = obj_list[[i]] %>% subset(condition == "Control_KO" | condition == "Control_WT")
  cell = cell_types[i]

  Idents(obj) <- "harmony_clusters"
  DEG_df = FindAllMarkers(obj, only.pos = TRUE) 
  DEG_df %>% write.csv(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", Sys.Date(), "_", cell, "_FindAllMarkers.csv"))
  
  for(j in 1:length(unique(DEG_df$cluster))){
    clust = unique(DEG_df$cluster)[j]
    output_path <- file.path(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/"))
    enrichR_pathway(DEG_df %>% filter(cluster == clust), pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = output_path, title = paste0(cell, "_", clust))
  }
  }, error = function(e) {
    cat("An error occurred:", conditionMessage(e), "\n")
  })
}


#########################################################
# DEG analysis on Control_Ms4a4a_KO vs Control_Ms4a4a_WT
#########################################################

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
protein_coding_genes <- getBM(attributes = c('mgi_symbol'), 
                              filters = 'biotype', 
                              values = 'protein_coding', 
                              mart = ensembl)
protein_coding_genes <- protein_coding_genes$mgi_symbol
head(protein_coding_genes)
 
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

# Function to process a single cell type
process_DEG_by_cluster <- function(seurat_obj, celltype_name, protein_coding_genes, output_dir) {
  seurat_obj = seurat_obj %>% subset(condition == "Control_KO" | condition == "Control_WT")
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
    out_csv <- file.path(output_dir, paste0("3-DEG_analysis/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    
    # Write to CSV
    write.csv(markers, out_csv)
    
    #Volcano plot
    output_vol <- file.path(output_dir, paste0("4-Volcano_plot/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    alternative_volcano_plot(markers, "p_val_adj", "avg_log2FC", path = output_vol, title = paste0(celltype_name, "_KO_vs_WT_", clust), genotype = "KO vs WT", celltype = celltype_name)

    #EnrichR
    output_path <- file.path(output_dir, paste0("5-EnrichR/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
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



##################################################################################################################
##################################################################################################################
###                                                                                                           ####
###                                                Only female                                                ####
###                                                                                                           ####
##################################################################################################################
##################################################################################################################


#########################################################
# UMAP
#########################################################

obj_list = list(Micro, Astro, Excit, Inhib, Oligo, OPC, Other)
cell_types = c("0-Microglia", "1-Astrocytes", "2-Excitatory_neurons", "3-Inhibitory_neurons", "4-Oligodendrocyte", "5-OPC", "6-Other")

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/")

for(i in 1:length(obj_list)){
  obj = obj_list[[i]] %>% subset(sex == "F")
  cell = cell_types[i]
  png(paste0(cell, "/0-UMAP/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", cell, "_UMAP.png"), width = 6, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7) + 
    labs(title = paste0(cell, " Subclusters")) +
    theme(text = element_text(size=20, face="bold"))
  print(p)
  dev.off()
  png(paste0(cell, "/0-UMAP/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", cell, "_UMAP_KO_vs_WT.png"), width = 12, height = 5, units = "in", res = 300)
  p = DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters",label = TRUE, pt.size = 0.7, split.by = "condition") + 
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
  tryCatch({
    obj = obj_list[[i]] %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
    cell = cell_types[i]

    output_dir <- paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/1-Proportion/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/")
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
  }, error = function(e) {
    cat("An error occurred:", conditionMessage(e), "\n")
  })
}


# No significant cluster


#########################################################
# FindAllMarkers
#########################################################

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

for(i in 1:length(obj_list)){
  tryCatch({
  obj = obj_list[[i]] %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
  cell = cell_types[i]

  Idents(obj) <- "harmony_clusters"
  DEG_df = FindAllMarkers(obj, only.pos = TRUE) 
  DEG_df %>% write.csv(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", Sys.Date(), "_", cell, "_FindAllMarkers.csv"))
  
  for(j in 1:length(unique(DEG_df$cluster))){
    clust = unique(DEG_df$cluster)[j]
    output_path <- file.path(paste0("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/2-Subcluster_analysis/", cell, "/2-FindAllMarkers/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"))
    enrichR_pathway(DEG_df %>% filter(cluster == clust), pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = output_path, title = paste0(cell, "_", clust))
  }
  }, error = function(e) {
    cat("An error occurred:", conditionMessage(e), "\n")
  })
}


#########################################################
# DEG analysis on 5xFAD_Ms4a4a_KO vs 5xFAD_Ms4a4a_WT
#########################################################

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
protein_coding_genes <- getBM(attributes = c('mgi_symbol'), 
                              filters = 'biotype', 
                              values = 'protein_coding', 
                              mart = ensembl)
protein_coding_genes <- protein_coding_genes$mgi_symbol
head(protein_coding_genes)
 
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")

# Function to process a single cell type
process_DEG_by_cluster <- function(seurat_obj, celltype_name, protein_coding_genes, output_dir) {
  seurat_obj = seurat_obj %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
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
    out_csv <- file.path(output_dir, paste0("3-DEG_analysis/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    
    # Write to CSV
    write.csv(markers, out_csv)
    
    #Volcano plot
    output_vol <- file.path(output_dir, paste0("4-Volcano_plot/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
    alternative_volcano_plot(markers, "p_val_adj", "avg_log2FC", path = output_vol, title = paste0(celltype_name, "_KO_vs_WT_", clust), genotype = "KO vs WT", celltype = celltype_name)

    #EnrichR
    output_path <- file.path(output_dir, paste0("5-EnrichR/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/", format(Sys.Date(), "%Y-%m-%d"), " KO_vs_WT_", celltype_name, "_", clust, ".csv"))
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

