conda activate seurat_v5
R

library(qs)
library(Seurat)
library(dplyr)
library(ggrepel)
library(stringr)
library(Azimuth)
library(scales)

merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj$fAD = ifelse(str_detect(merged_obj$orig.ident, "Exp12"), "Control", "5xFAD")
merged_obj$condition = paste0(merged_obj$fAD, "_", merged_obj$genotype)
merged_obj@meta.data %>% head
merged_obj %>% qsave("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")


##################################################################################################################
##################################################################################################################
###                                                                                                           ####
###                                       All the genotype and sex                                            ####
###                                                                                                           ####
##################################################################################################################
##################################################################################################################


#########################################################
# Calculate Cell type Proportion by Ms4a4a KO
#########################################################


df = merged_obj@meta.data
df %>% group_by(condition) %>% summarise(n = n_distinct(orig.ident))

#   condition      n
#   <chr>      <int>
# 1 5xFAD_KO       6
# 2 5xFAD_WT       6
# 3 Control_KO     2
# 4 Control_WT     2

df = df %>% group_by(orig.ident, manual_annotation, condition) %>% summarise(n = n()) %>% ungroup
#For each animal, calculate the percentage of each cell type
df = df %>% group_by(orig.ident) %>% mutate(percentage = n/sum(n)) %>% ungroup
df %>% write.csv("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/0-Proportion/2025-4-13 merged_obj_proportion.csv")

for(cell in unique(df$manual_annotation)){
  print(cell)
  test = df %>% filter(manual_annotation == cell)
  res = pairwise.t.test(test$percentage, test$condition, paired = F, p.adjust.method = "none", pool.sd = F)
  print(res$p.value)
}

# [1] "Astro"
#              5xFAD_KO    5xFAD_WT Control_KO
# 5xFAD_WT   0.28803624          NA         NA
# Control_KO 0.11071216 0.103259050         NA
# Control_WT 0.03786349 0.002116631  0.2633535
# [1] "Excit"
#              5xFAD_KO   5xFAD_WT Control_KO
# 5xFAD_WT   0.88497314         NA         NA
# Control_KO 0.10876095 0.12561506         NA
# Control_WT 0.01331147 0.01602496  0.2984783
# [1] "Inhib"
#              5xFAD_KO    5xFAD_WT Control_KO
# 5xFAD_WT   0.06175619          NA         NA
# Control_KO 0.00989233 0.001589737         NA
# Control_WT 0.18302584 0.026815214  0.4588053
# [1] "Micro"
#                5xFAD_KO     5xFAD_WT Control_KO
# 5xFAD_WT   0.5708385318           NA         NA
# Control_KO 0.0003720972 0.0002592610         NA
# Control_WT 0.0002884459 0.0001839507  0.6569011
# [1] "OPC"
#              5xFAD_KO   5xFAD_WT Control_KO
# 5xFAD_WT   0.69580927         NA         NA
# Control_KO 0.07674667 0.04597657         NA
# Control_WT 0.12813452 0.07297945   0.131758
# [1] "Oligo"
#             5xFAD_KO   5xFAD_WT Control_KO
# 5xFAD_WT   0.2953395         NA         NA
# Control_KO 0.4458005 0.01181260         NA
# Control_WT 0.4245877 0.01176217  0.8710796
# [1] "Other"
#              5xFAD_KO   5xFAD_WT Control_KO
# 5xFAD_WT   0.84003753         NA         NA
# Control_KO 0.10595640 0.07982427         NA
# Control_WT 0.03017444 0.01505198  0.1708488
# [1] "Unknown"
#              5xFAD_KO    5xFAD_WT Control_KO
# 5xFAD_WT   0.21675884          NA         NA
# Control_KO 0.00219224 0.004194984         NA
# Control_WT 0.01749892 0.005211480  0.7980115

df$manual_annotation <- factor(df$manual_annotation,
                               levels = unique(df$manual_annotation))
df$condition <- factor(df$condition,
                       levels = c("Control_WT", "Control_KO", "5xFAD_WT", "5xFAD_KO"))

# Plot
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/0-Proportion/2025-4-13 merged_obj_proportion.png", width = 800, height = 400)
ggplot(df, aes(x = manual_annotation, y = percentage, fill = condition, color = condition)) +
  geom_bar(stat = "summary", fun = "mean", 
           position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  geom_point(aes(group = condition), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.8, alpha = 0.9) +
  labs(x = "Cell Type", y = "Relative Abundance (Mean %)", fill = "Condition", color = "Condition", title = "All animal (Control/5xFAD, Ms4a4a WT/KO)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
        
dev.off()

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

#---------------------------------------------------------------- 5xFAD_Ms4a4a KO vs 5xFAD_Ms4a4a WT ----------------------------------------------------------------

merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj = merged_obj %>% subset(condition == "5xFAD_KO" | condition == "5xFAD_WT")
Idents(merged_obj) <- "condition"

Ext.de <- merged_obj %>% subset(manual_annotation == "Excit") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Inh.de <- merged_obj %>% subset(manual_annotation == "Inhib") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Ast.de <- merged_obj %>% subset(manual_annotation == "Astro") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Mic.de <- merged_obj %>% subset(manual_annotation == "Micro") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Opc.de <- merged_obj %>% subset(manual_annotation == "OPC") %>%   FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE)   
Oli.de <- merged_obj %>% subset(manual_annotation == "Oligo") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Oth.de <- merged_obj %>% subset(manual_annotation == "Other") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 

Ext.de = Ext.de %>% mutate(gene_biotype = ifelse(row.names(Ext.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ext.de))
Inh.de = Inh.de %>% mutate(gene_biotype = ifelse(row.names(Inh.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Inh.de))
Ast.de = Ast.de %>% mutate(gene_biotype = ifelse(row.names(Ast.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ast.de))
Mic.de = Mic.de %>% mutate(gene_biotype = ifelse(row.names(Mic.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.de))
Opc.de = Opc.de %>% mutate(gene_biotype = ifelse(row.names(Opc.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Opc.de))
Oli.de = Oli.de %>% mutate(gene_biotype = ifelse(row.names(Oli.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oli.de))
Oth.de = Oth.de %>% mutate(gene_biotype = ifelse(row.names(Oth.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oth.de))

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/2-DEG_analysis/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/")
Ext.de %>% write.csv("2025-4-13 KO_vs_WT_Ext.csv")
Inh.de %>% write.csv("2025-4-13 KO_vs_WT_Inh.csv")
Ast.de %>% write.csv("2025-4-13 KO_vs_WT_Ast.csv")
Mic.de %>% write.csv("2025-4-13 KO_vs_WT_Mic.csv")
Oli.de %>% write.csv("2025-4-13 KO_vs_WT_Oli.csv")
Opc.de %>% write.csv("2025-4-13 KO_vs_WT_Opc.csv")
Oth.de %>% write.csv("2025-4-13 KO_vs_WT_Oth.csv")

Ext.de = read.csv("2025-4-13 KO_vs_WT_Ext.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Inh.de = read.csv("2025-4-13 KO_vs_WT_Inh.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Ast.de = read.csv("2025-4-13 KO_vs_WT_Ast.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.de = read.csv("2025-4-13 KO_vs_WT_Mic.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oli.de = read.csv("2025-4-13 KO_vs_WT_Oli.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Opc.de = read.csv("2025-4-13 KO_vs_WT_Opc.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oth.de = read.csv("2025-4-13 KO_vs_WT_Oth.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/3-Volcano_plot/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"

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
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/4-Pathway/All_cell_type/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"

ExtPath = enrichR_pathway(Ext.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Excitatory_neuron_KO_vs_WT")
InhPath = enrichR_pathway(Inh.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Inhibitory_neuron_KO_vs_WT")
AstPath = enrichR_pathway(Ast.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Astrocyte_KO_vs_WT")
MicPath = enrichR_pathway(Mic.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Microglia_KO_vs_WT")
OliPath = enrichR_pathway(Oli.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Oligodendrocyte_KO_vs_WT")
OpcPath = enrichR_pathway(Opc.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "OPC_KO_vs_WT")
OthPath = enrichR_pathway(Oth.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Other_KO_vs_WT")

#---------------------------------------------------------------- Control_Ms4a4a KO vs Control_Ms4a4a WT ----------------------------------------------------------------

merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj = merged_obj %>% subset(condition == "Control_KO" | condition == "Control_WT")
Idents(merged_obj) <- "condition"

Ext.de <- merged_obj %>% subset(manual_annotation == "Excit") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 
Inh.de <- merged_obj %>% subset(manual_annotation == "Inhib") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 
Ast.de <- merged_obj %>% subset(manual_annotation == "Astro") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 
Mic.de <- merged_obj %>% subset(manual_annotation == "Micro") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 
Opc.de <- merged_obj %>% subset(manual_annotation == "OPC") %>%   FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE)   
Oli.de <- merged_obj %>% subset(manual_annotation == "Oligo") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 
Oth.de <- merged_obj %>% subset(manual_annotation == "Other") %>% FindMarkers(ident.1 = "Control_KO", ident.2 = "Control_WT", verbose = FALSE) 

Ext.de = Ext.de %>% mutate(gene_biotype = ifelse(row.names(Ext.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ext.de))
Inh.de = Inh.de %>% mutate(gene_biotype = ifelse(row.names(Inh.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Inh.de))
Ast.de = Ast.de %>% mutate(gene_biotype = ifelse(row.names(Ast.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ast.de))
Mic.de = Mic.de %>% mutate(gene_biotype = ifelse(row.names(Mic.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.de))
Opc.de = Opc.de %>% mutate(gene_biotype = ifelse(row.names(Opc.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Opc.de))
Oli.de = Oli.de %>% mutate(gene_biotype = ifelse(row.names(Oli.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oli.de))
Oth.de = Oth.de %>% mutate(gene_biotype = ifelse(row.names(Oth.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oth.de))

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/2-DEG_analysis/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT")
Ext.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Ext.csv")
Inh.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Inh.csv")
Ast.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Ast.csv")
Mic.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Mic.csv")
Oli.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Oli.csv")
Opc.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Opc.csv")
Oth.de %>% write.csv("2025-4-15_Control_KO_vs_Control_WT_Oth.csv")

Ext.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Ext.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Inh.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Inh.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Ast.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Ast.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Mic.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oli.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Oli.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Opc.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Opc.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oth.de = read.csv("2025-4-15_Control_KO_vs_Control_WT_Oth.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/3-Volcano_plot/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/"

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
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/4-Pathway/All_cell_type/Control_Ms4a4a_KO_vs_Control_Ms4a4a_WT/"

ExtPath = enrichR_pathway(Ext.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Excitatory_neuron_KO_vs_WT")
InhPath = enrichR_pathway(Inh.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Inhibitory_neuron_KO_vs_WT")
AstPath = enrichR_pathway(Ast.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Astrocyte_KO_vs_WT")
MicPath = enrichR_pathway(Mic.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Microglia_KO_vs_WT")
OliPath = enrichR_pathway(Oli.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Oligodendrocyte_KO_vs_WT")
OpcPath = enrichR_pathway(Opc.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "OPC_KO_vs_WT")
OthPath = enrichR_pathway(Oth.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Other_KO_vs_WT")

##################################################################################################################
##################################################################################################################
###                                                                                                           ####
###                                                Only female                                                ####
###                                                                                                           ####
##################################################################################################################
##################################################################################################################

#########################################################
# Calculate Cell type Proportion by Ms4a4a KO
#########################################################

df = merged_obj@meta.data %>% filter(sex == "F")
df %>% group_by(condition) %>% summarise(n = n_distinct(orig.ident))

#   condition     n
#   <chr>     <int>
# 1 5xFAD_KO      3
# 2 5xFAD_WT      3

df = df %>% group_by(orig.ident, manual_annotation, condition) %>% summarise(n = n()) %>% ungroup
#For each animal, calculate the percentage of each cell type
df = df %>% group_by(orig.ident) %>% mutate(percentage = n/sum(n)) %>% ungroup
df %>% write.csv("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/0-Proportion/2025-4-13 merged_obj_proportion_only_female.csv")

for(cell in unique(df$manual_annotation)){
  print(cell)
  test = df %>% filter(manual_annotation == cell)
  res = pairwise.t.test(test$percentage, test$condition, paired = F, p.adjust.method = "none", pool.sd = F)
  print(res$p.value)
 }

# [1] "Astro"
#           5xFAD_KO
# 5xFAD_WT 0.3811053
# [1] "Excit"
#          5xFAD_KO
# 5xFAD_WT 0.727439
# [1] "Inhib"
#           5xFAD_KO
# 5xFAD_WT 0.2612412
# [1] "Micro"
#           5xFAD_KO
# 5xFAD_WT 0.5000948
# [1] "OPC"
#           5xFAD_KO
# 5xFAD_WT 0.9619964
# [1] "Oligo"
#           5xFAD_KO
# 5xFAD_WT 0.8632216
# [1] "Other"
#           5xFAD_KO
# 5xFAD_WT 0.6507127
# [1] "Unknown"
#           5xFAD_KO
# 5xFAD_WT 0.4985041

df$manual_annotation <- factor(df$manual_annotation,
                               levels = unique(df$manual_annotation))
df$condition <- factor(df$condition,
                       levels = c("Control_WT", "Control_KO", "5xFAD_WT", "5xFAD_KO"))

# Plot
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/0-Proportion/2025-4-13 merged_obj_proportion_only_female.png", width = 800, height = 400)
ggplot(df, aes(x = manual_annotation, y = percentage, fill = condition, color = condition)) +
  geom_bar(stat = "summary", fun = "mean", 
           position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  geom_point(aes(group = condition), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.8, alpha = 0.9) +
  labs(x = "Cell Type", y = "Relative Abundance (Mean %)", fill = "Condition", color = "Condition", title = "All animal (Control/5xFAD, Ms4a4a WT/KO)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
        
dev.off()


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

merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj = merged_obj %>% subset(sex == "F")
Idents(merged_obj) <- "condition"

Ext.de <- merged_obj %>% subset(manual_annotation == "Excit") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Inh.de <- merged_obj %>% subset(manual_annotation == "Inhib") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Ast.de <- merged_obj %>% subset(manual_annotation == "Astro") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Mic.de <- merged_obj %>% subset(manual_annotation == "Micro") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Opc.de <- merged_obj %>% subset(manual_annotation == "OPC") %>%   FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE)   
Oli.de <- merged_obj %>% subset(manual_annotation == "Oligo") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 
Oth.de <- merged_obj %>% subset(manual_annotation == "Other") %>% FindMarkers(ident.1 = "5xFAD_KO", ident.2 = "5xFAD_WT", verbose = FALSE) 

Ext.de = Ext.de %>% mutate(gene_biotype = ifelse(row.names(Ext.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ext.de))
Inh.de = Inh.de %>% mutate(gene_biotype = ifelse(row.names(Inh.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Inh.de))
Ast.de = Ast.de %>% mutate(gene_biotype = ifelse(row.names(Ast.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Ast.de))
Mic.de = Mic.de %>% mutate(gene_biotype = ifelse(row.names(Mic.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Mic.de))
Opc.de = Opc.de %>% mutate(gene_biotype = ifelse(row.names(Opc.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Opc.de))
Oli.de = Oli.de %>% mutate(gene_biotype = ifelse(row.names(Oli.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oli.de))
Oth.de = Oth.de %>% mutate(gene_biotype = ifelse(row.names(Oth.de) %in% protein_coding_genes, "protein_coding_genes", NA), gene = rownames(Oth.de))

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/2-DEG_analysis/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT")
Ext.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Ext.csv")
Inh.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Inh.csv")
Ast.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Ast.csv")
Mic.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Mic.csv")
Oli.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Oli.csv")
Opc.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Opc.csv")
Oth.de %>% write.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Oth.csv")

Ext.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Ext.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Inh.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Inh.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Ast.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Ast.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Mic.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Mic.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oli.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Oli.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Opc.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Opc.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")
Oth.de = read.csv("2025-4-15_5xFAD_KO_vs_5xFAD_WT_Oth.csv", row.names = 1) %>% filter(gene_biotype == "protein_coding_genes")

source("/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Tool_code/01_Volcano_plot.r")
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/3-Volcano_plot/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"

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
Path = "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/1-Major_celltype_analysis/4-Pathway/Female_only/5xFAD_Ms4a4a_KO_vs_5xFAD_Ms4a4a_WT/"

ExtPath = enrichR_pathway(Ext.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Excitatory_neuron_KO_vs_WT")
InhPath = enrichR_pathway(Inh.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Inhibitory_neuron_KO_vs_WT")
AstPath = enrichR_pathway(Ast.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Astrocyte_KO_vs_WT")
MicPath = enrichR_pathway(Mic.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Microglia_KO_vs_WT")
OliPath = enrichR_pathway(Oli.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Oligodendrocyte_KO_vs_WT")
OpcPath = enrichR_pathway(Opc.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "OPC_KO_vs_WT")
OthPath = enrichR_pathway(Oth.de, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = "gene", path = Path, title = "Other_KO_vs_WT")

