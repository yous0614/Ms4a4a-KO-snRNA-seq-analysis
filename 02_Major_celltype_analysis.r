conda activate seurat_v5
R

library(qs)
library(Seurat)
library(dplyr)
library(ggrepel)
library(stringr)
library(Azimuth)

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


