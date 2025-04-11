conda activate seurat_v5
R

library(qs)
library(Seurat)
library(dplyr)
library(ggrepel)
library(stringr)
library(Azimuth)
library(scales)

#########################################################
# Calculate Cell type Proportion by Ms4a4a KO
#########################################################

merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
df = merged_obj@meta.data
df %>% group_by(genotype) %>% summarise(n = n_distinct(orig.ident))

# 1 WT           8
# 2 KO           8

df = df %>% group_by(orig.ident, manual_annotation, genotype) %>% summarise(n = n()) %>% ungroup
#For each animal, calculate the percentage of each cell type
df = df %>% group_by(orig.ident) %>% mutate(percentage = n/sum(n)) %>% ungroup
df %>% write.csv("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/0-Proportion/2025-4-6 merged_obj_proportion.csv")

for(cell in unique(df$manual_annotation)){
  print(cell)
  test = df %>% filter(manual_annotation == cell)
  res = pairwise.t.test(test$percentage, test$genotype, paired = F, p.adjust.method = "none", pool.sd = F)
  print(res$p.value)
}

# [1] "Astro"
#           WT
# KO 0.2272964
# [1] "Excit"
#          WT
# KO 0.664011
# [1] "Inhib"
#           WT
# KO 0.1150971
# [1] "Micro"
#           WT
# KO 0.7658652
# [1] "OPC"
#           WT
# KO 0.6583632
# [1] "Oligo"
#           WT
# KO 0.3188048
# [1] "Other"
#           WT
# KO 0.6923117
# [1] "Unknown"
#           WT
# KO 0.4082052

#########################################################
# FindAllMarkers for all the major celltypes
#########################################################

library(qs)
library(dplyr)

Idents(merged_obj) = "manual_annotation"
FindAllMarkers(merged_obj, only.pos = T) %>% write.csv("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Analysis/1-FindFallMarkers/2025-4-6 merged_obj_FindAllMarkers.csv")


#########################################################
# DEG analysis and Volcano plot for major cell types
#########################################################

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
protein_coding_genes <- getBM(attributes = c('mgi_symbol'), 
                              filters = 'biotype', 
                              values = 'protein_coding', 
                              mart = ensembl)
protein_coding_genes <- protein_coding_genes$mgi_symbol
head(protein_coding_genes)

merged_obj = merged_obj 
Idents(merged) <- "genotype"

Ext.de <- merged_obj_coding %>% subset(manual_annotation == "Excit") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Ext.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ext.de))
Inh.de <- merged_obj_coding %>% subset(manual_annotation == "Inhib") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Inh.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Inh.de))
Ast.de <- merged_obj_coding %>% subset(manual_annotation == "Astro") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Ast.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ast.de))
Mic.de <- merged_obj_coding %>% subset(manual_annotation == "Micro") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Mic.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.de))
Opc.de <- merged_obj_coding %>% subset(manual_annotation == "OPC") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE)   %>% mutate(gene_biotype = ifelse(row.names(Opc.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Opc.de))
Oli.de <- merged_obj_coding %>% subset(manual_annotation == "Oligo") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Oli.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oli.de))
Oth.de <- merged_obj_coding %>% subset(manual_annotation == "Other") %>% FindMarkers(ident.1 = "KO", ident.2 = "WT", verbose = FALSE) %>% mutate(gene_biotype = ifelse(row.names(Oth.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oth.de))

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/2-DEG_analysis")
Ext.de %>% write.csv("2025-4-6 KO_vs_WT_Ext.csv")
Inh.de %>% write.csv("2025-4-6 KO_vs_WT_Inh.csv")
Ast.de %>% write.csv("2025-4-6 KO_vs_WT_Ast.csv")
Mic.de %>% write.csv("2025-4-6 KO_vs_WT_Mic.csv")
Oli.de %>% write.csv("2025-4-6 KO_vs_WT_Oli.csv")
Opc.de %>% write.csv("2025-4-6 KO_vs_WT_Opc.csv")
Oth.de %>% write.csv("2025-4-6 KO_vs_WT_Oth.csv")

Ext.de = read.csv("2025-4-6 KO_vs_WT_Ext.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Inh.de = read.csv("2025-4-6 KO_vs_WT_Inh.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Ast.de = read.csv("2025-4-6 KO_vs_WT_Ast.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.de = read.csv("2025-4-6 KO_vs_WT_Mic.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oli.de = read.csv("2025-4-6 KO_vs_WT_Oli.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Opc.de = read.csv("2025-4-6 KO_vs_WT_Opc.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oth.de = read.csv("2025-4-6 KO_vs_WT_Oth.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")

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
# Pathway analysis for major cell types
#########################################################

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/02_EnrichR.r")
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/4-Pathway/"

ExtPath = enrichR_pathway(Ext.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Excitatory_neuron_KO_vs_WT")
InhPath = enrichR_pathway(Inh.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Inhibitory_neuron_KO_vs_WT")
AstPath = enrichR_pathway(Ast.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Astrocyte_KO_vs_WT")
MicPath = enrichR_pathway(Mic.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Microglia_KO_vs_WT")
OliPath = enrichR_pathway(Oli.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Oligodendrocyte_KO_vs_WT")
OpcPath = enrichR_pathway(Opc.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "OPC_KO_vs_WT")
OthPath = enrichR_pathway(Oth.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Other_KO_vs_WT")


