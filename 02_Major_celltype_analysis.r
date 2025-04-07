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
DF = Trem2PIC_microglia@meta.data
DF <- DF %>%
  mutate(microglia_subcluster = case_when(
    harmony_clusters == 0 ~ "Mic.0",
    harmony_clusters == 1 ~ "Mic.1",
    harmony_clusters == 2 ~ "Mic.2",
    harmony_clusters == 3 ~ "BAM",
    harmony_clusters == 4 ~ "Mic.3",
    TRUE ~ as.character(harmony_clusters)
  ))
Trem2PIC_microglia = AddMetaData(Trem2PIC_microglia, DF$microglia_subcluster, col.name = "microglia_subcluster")

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

#########################################################
#Calculate Proportional significancy
#########################################################

#Only M, WT/KO
df = Trem2PIC_microglia@meta.data %>% filter(Sex == "M") %>% filter(Genotype %in% c("WT", "KO"))
df$Condition = paste(df$Genotype, df$Treatment, sep = "_")
df %>% group_by(Condition) %>% summarise(n = n_distinct(orig.ident))
# 1 KO_POLY_I-C     3
# 2 KO_VEH          2
# 3 WT_POLY_I-C     4
# 4 WT_VEH          2
df = df %>% group_by(orig.ident, microglia_subcluster, Condition) %>% summarise(n = n()) %>% ungroup
#For each animal, calculate the percentage of each cell type
df = df %>% group_by(orig.ident) %>% mutate(percentage = n/sum(n)) %>% ungroup
df %>% write.csv("/home/yous/TREM2XSCZ/2025-3-1 Do it all over again/2-result/0-Process/Proportion/2025-3-1 Trem2PIC_microglia_proportion.csv")

for(cell in unique(df$microglia_subcluster)){
  print(cell)
  test = df %>% filter(microglia_subcluster == cell)
  res = pairwise.t.test(test$percentage, test$Condition, paired = F, p.adjust.method = "none", pool.sd = F)
  print(res)
}

# [1] "Mic.0"

#         Pairwise comparisons using t tests with non-pooled SD 

# data:  test$percentage and test$Condition 

#             KO_POLY_I-C KO_VEH WT_POLY_I-C
# KO_VEH      0.60        -      -          
# WT_POLY_I-C 0.67        0.81   -          
# WT_VEH      0.39        0.42   0.42       

# P value adjustment method: none 
# [1] "Mic.1"

#         Pairwise comparisons using t tests with non-pooled SD 

# data:  test$percentage and test$Condition 

#             KO_POLY_I-C KO_VEH WT_POLY_I-C
# KO_VEH      0.55        -      -          
# WT_POLY_I-C 0.89        0.61   -          
# WT_VEH      0.41        0.43   0.41       

# P value adjustment method: none 
# [1] "Mic.2"

#         Pairwise comparisons using t tests with non-pooled SD 

# data:  test$percentage and test$Condition 

#             KO_POLY_I-C KO_VEH WT_POLY_I-C
# KO_VEH      0.746       -      -          
# WT_POLY_I-C 0.076       0.054  -          
# WT_VEH      0.272       0.241  0.960      

# P value adjustment method: none 
# [1] "Mic.3"

#         Pairwise comparisons using t tests with non-pooled SD 

# data:  test$percentage and test$Condition 

#             KO_POLY_I-C KO_VEH WT_POLY_I-C
# KO_VEH      0.91        -      -          
# WT_POLY_I-C 0.83        0.68   -          
# WT_VEH      0.94        0.96   0.67       

# P value adjustment method: none 
# [1] "BAM"

#         Pairwise comparisons using t tests with non-pooled SD 

# data:  test$percentage and test$Condition 

#             KO_POLY_I-C KO_VEH WT_POLY_I-C
# KO_VEH      0.45        -      -          
# WT_POLY_I-C 0.68        0.09   -          
# WT_VEH      0.90        0.35   0.45       

# P value adjustment method: none 