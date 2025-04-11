conda activate seurat_v5
R

library(qs)
library(Seurat)
library(dplyr)
library(ggrepel)
library(stringr)
library(Azimuth)

#########################################################
# Make Seurat Object
#########################################################

today = Sys.Date()

setwd("/home/data/WashU_Data/Multiome_Data/2023_12_08_Karch_SR003027/")
file_path = list.dirs(".", recursive = TRUE) %>% grep("filtered_feature_bc_matrix", ., value = TRUE)
file_name = gsub("/outs/filtered_feature_bc_matrix", "", file_path) %>% gsub("./SR003027_", "", .)
datasets = list()

for(i in 1:length(file_path)){
print(file_path[i])
obj = Read10X(file_path[i])
obj = CreateSeuratObject(counts = obj$`Gene Expression`, project = file_name[i])
obj = RunAzimuth(obj, reference = "/home/yous/github/Ms4a4a-KO-snRNA-seq-analysis/Reference/Azimuth-Mice_Motor_Cortex/") 
obj = RenameCells(obj, add.cell.id = paste0("Sample_", file_name[i]))
datasets[i] = obj
}

length(datasets)
merged_obj = merge(x = datasets[[1]], y = datasets[2:16])

# Adding Ms4a4a WT/KO and sex
df = merged_obj@meta.data
df$orig.ident %>% table
df$genotype = ifelse(str_detect(df$orig.ident, "54|57|75|53|35|70|WT"), "WT", 
                     ifelse(str_detect(df$orig.ident, "60|47|73|42|59|66|KO"), "KO", NA))
df$sex = ifelse(str_detect(df$orig.ident, "54|57|75|60|47|73|WT|KO"), "M", "F")
df$sex = factor(df$sex, levels = c("M", "F"))
df$genotype = factor(df$genotype, levels = c("WT", "KO"))
merged_obj = AddMetaData(merged_obj, metadata = df)

# Processing the merged object
merged_obj = merged_obj %>%
  PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^Rps|^Rpl", col.name = "percent.ribo")

Idents(merged_obj) = "orig.ident"

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/0-Filter_low_qualities_cells/")
png(paste0(today, "_MS4A_KO_mice_before_QC_merged.png"), width = 15, height = 8, units = "in", res = 300)
VlnPlot(merged_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 4)
dev.off()

merged_obj = merged_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 0.2) %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#####################################################################
# Filter QC by lower 1% UMI counts per donor and <5% mito 
#####################################################################

for(i in unique(merged_obj$orig.ident)){
    print(i)
    umi_counts = merged_obj %>% subset(orig.ident == i) %>% .$nCount_RNA %>% log10()
    
    q = quantile(umi_counts, probs = c(0.01, 0.02, 0.03))
    q_labels = c("1%", "2%", "3%")
    q_colors = c("red", "green", "blue")
    
    # Make histogram
    setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/0-Filter_low_qualities_cells/")
    png(paste0("UMI_histogram_", i, ".png"), width = 10, height = 6, units = "in", res = 300)
    p = ggplot(data.frame(UMI = umi_counts), aes(x = UMI)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
    labs(x = "UMI count per cell", y = "Number of cells", title = paste0(i)) +
    theme_minimal() +
    geom_vline(xintercept = q[1], color = q_colors[1], linetype = "dashed", size = 1) +
    geom_vline(xintercept = q[2], color = q_colors[2], linetype = "dashed", size = 1) +
    geom_vline(xintercept = q[3], color = q_colors[3], linetype = "dashed", size = 1) +
    annotate("text", x = q, y = Inf, label = q_labels, color = q_colors,
            vjust = -0.5, hjust = -0.1, angle = 90, size = 3)
    print(p)
    dev.off()
}

merged_obj$QC_Filtering <- NA  # initialize column

for(i in unique(merged_obj$orig.ident)) {
    print(i)
    
    # Subset the cells from the current sample
    cells_i <- WhichCells(merged_obj %>% subset(orig.ident = i))  # safer indexing
    umi_counts <- merged_obj$nCount_RNA[cells_i]

    # Calculate low quantile thresholds (e.g., bottom 1%)
    q <- quantile(umi_counts, probs = c(0.01))

    # Set filtering status
    merged_obj$QC_Filtering[cells_i] <- ifelse(umi_counts < q, "remove", "retain")
}


merged_obj$QC_Filtering %>% table
# remove retain 
#   1571 155767 

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/0-Filter_low_qualities_cells")
png(paste0(today, "_UMAP_MS4A_KO_mice_QC_filtering.png"), width = 15, height = 6, units = "in", res = 300)
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("predicted.subclass", "QC_Filtering"), label = TRUE) 
dev.off()
png(paste0(today, "_UMAP_MS4A_KO_mice_QC_filtering_per_donor.png"), width = 15, height = 15, units = "in", res = 300)
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "QC_Filtering", split.by = "orig.ident", ncol = 4, pt.size = 2)
dev.off()

qsave(merged_obj, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/01_MS4A_KO_mice_unintegrated.qs")

# Before filtering
merged_obj
merged_obj@meta.data$nCount_RNA %>% summary
merged_obj@meta.data$nFeature_RNA %>% summary
merged_obj@meta.data$percent.mt %>% summary
merged_obj@meta.data$percent.ribo %>% summary
merged_obj@meta.data$predicted.subclass %>% table

# 32473 features across 157338 samples within 5 assays 
# Active assay: RNA (32285 features, 2000 variable features)
#  33 layers present: counts.5XFAD4A-35, counts.5XFAD4A-42, counts.5XFAD4A-47, counts.5XFAD4A-53, counts.5XFAD4A-54, counts.5XFAD4A-57, counts.5XFAD4A-59, counts.5XFAD4A-60, counts.5XFAD4A-66, counts.5XFAD4A-70, counts.5XFAD4A-73, counts.5XFAD4A-75, counts.Exp12_KO1_Hipp, counts.Exp12_KO2_Hipp, counts.Exp12_WT1_Hipp, counts.Exp12_WT2_Hipp, data.5XFAD4A-35, data.5XFAD4A-42, data.5XFAD4A-47, data.5XFAD4A-53, data.5XFAD4A-54, data.5XFAD4A-57, data.5XFAD4A-59, data.5XFAD4A-60, data.5XFAD4A-66, data.5XFAD4A-70, data.5XFAD4A-73, data.5XFAD4A-75, data.Exp12_KO1_Hipp, data.Exp12_KO2_Hipp, data.Exp12_WT1_Hipp, data.Exp12_WT2_Hipp, scale.data
#  4 other assays present: prediction.score.class, prediction.score.cluster, prediction.score.subclass, prediction.score.cross_species_cluster
#  2 dimensional reductions calculated: pca, umap.unintegrated
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       4    1262    3545    6212    8469  704804 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       4     825    1879    2196    3227   15718 
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.00000  0.01071  0.06882  0.17561  0.18975 13.70968 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.2776  0.4170  0.5856  0.6662 12.3426 
# .
#      Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT 
#       9949        123      48411      13684       9288       2323       2549 
#      L6 IT L6 IT Car3        L6b      Lamp5      Meis2  Micro-PVM      Oligo 
#       3582        791       1541       3089        452      17959      28265 
#        OPC       Peri      Pvalb       Sncg        Sst  Sst Chodl        Vip 
#       5142        206       2685       1216       2826        141       1709 
#       VLMC 
#       1407 

merged_obj = merged_obj %>% subset(QC_Filtering == "retain" & percent.mt < 5)

# After filtering
merged_obj
merged_obj@meta.data$nCount_RNA %>% summary
merged_obj@meta.data$nFeature_RNA %>% summary
merged_obj@meta.data$percent.mt %>% summary
merged_obj@meta.data$percent.ribo %>% summary
merged_obj@meta.data$predicted.subclass %>% table

# 32473 features across 155679 samples within 5 assays 
# Active assay: RNA (32285 features, 2000 variable features)
#  33 layers present: counts.5XFAD4A-35, counts.5XFAD4A-42, counts.5XFAD4A-47, counts.5XFAD4A-53, counts.5XFAD4A-54, counts.5XFAD4A-57, counts.5XFAD4A-59, counts.5XFAD4A-60, counts.5XFAD4A-66, counts.5XFAD4A-70, counts.5XFAD4A-73, counts.5XFAD4A-75, counts.Exp12_KO1_Hipp, counts.Exp12_KO2_Hipp, counts.Exp12_WT1_Hipp, counts.Exp12_WT2_Hipp, data.5XFAD4A-35, data.5XFAD4A-42, data.5XFAD4A-47, data.5XFAD4A-53, data.5XFAD4A-54, data.5XFAD4A-57, data.5XFAD4A-59, data.5XFAD4A-60, data.5XFAD4A-66, data.5XFAD4A-70, data.5XFAD4A-73, data.5XFAD4A-75, data.Exp12_KO1_Hipp, data.Exp12_KO2_Hipp, data.Exp12_WT1_Hipp, data.Exp12_WT2_Hipp, scale.data
#  4 other assays present: prediction.score.class, prediction.score.cluster, prediction.score.subclass, prediction.score.cross_species_cluster
#  2 dimensional reductions calculated: pca, umap.unintegrated
# > merged_obj@meta.data$nCount_RNA %>% summary
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     246    1287    3650    6275    8537  704804 
# > merged_obj@meta.data$nFeature_RNA %>% summary
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     150     840    1919    2217    3241   15718 
# > merged_obj@meta.data$percent.mt %>% summary
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01163 0.06949 0.17092 0.18927 4.96894 
# > merged_obj@meta.data$percent.ribo %>% summary
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.2779  0.4159  0.5808  0.6614 12.3426 
# > merged_obj@meta.data$predicted.subclass %>% table
# .
#      Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT 
#       9936        108      48346      13684       9154       2323       2548 
#      L6 IT L6 IT Car3        L6b      Lamp5      Meis2  Micro-PVM      Oligo 
#       3196        791       1541       3084        451      17450      27771 
#        OPC       Peri      Pvalb       Sncg        Sst  Sst Chodl        Vip 
#       5142        193       2685       1216       2826        141       1708 
#       VLMC 
#       1385 

merged_obj = JoinLayers(merged_obj)
qsave(merged_obj, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/02_MS4A_KO_mice_filter_QC.qs")


#########################################################
# Doublet removal
#########################################################

library(scDblFinder)
library(BiocParallel)

sce = as.SingleCellExperiment(merged_obj)
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3))
doublet_class <- sce$scDblFinder.class
merged_obj$scDblFinder.class <- doublet_class

table(sce$scDblFinder.class)
# singlet doublet 
#  139516   16163 

library(DoubletFinder)
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep(merged_obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
pdf("output.pdf")  # Redirects graphical output to a PDF file
bcmvn <- find.pK(sweep.stats)
dev.off()  # Closes the graphics device

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(merged_obj@meta.data$ClusteringResults)           ## ex: annotations <- merged_obj@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(merged_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
merged_obj <- doubletFinder(merged_obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
merged_obj@meta.data %>% head
merged_obj@meta.data$DF.classifications_0.25_0.09_11676 %>% table
# Doublet Singlet 
#   11676  144003

merged_obj@meta.data[,c("scDblFinder.class", "DF.classifications_0.25_0.09_11676")] %>% table

# scDblFinder.class Doublet Singlet
#           singlet    9807  129709
#           doublet    1869   14294


png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/MS4A_KO_mice_umap_scDblFinder.png", width = 15, height = 6, units = "in", res = 300)
merged_obj$scDblFinder.class = factor(merged_obj$scDblFinder.class, levels = rev(c("singlet", "doublet")))
DimPlot( merged_obj, reduction = "umap.unintegrated", group.by = c("predicted.subclass", "scDblFinder.class"), label = TRUE) +
  labs(title = "scDblFinder Classification") 
dev.off()

png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/MS4A_KO_mice_umap_DoubletFinder.png", width = 15, height = 6, units = "in", res = 300)
merged_obj$DF.classifications_0.25_0.09_11676 = factor(merged_obj$DF.classifications_0.25_0.09_11676, levels = rev(c("Singlet", "Doublet")))
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("predicted.subclass", "DF.classifications_0.25_0.09_11676"), label = TRUE) +
  labs(title = "DoubletFinder Classification") 
dev.off()

library(patchwork)
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/Microglia.png", width = 25, height = 10, units = "in", res = 300)
p1 = DimPlot(merged_obj %>% subset(predicted.subclass == "Micro-PVM"), reduction = "umap.unintegrated", label = TRUE, pt.size = 1.5)
# p2 = DimPlot(merged_obj %>% subset(predicted.subclass == "Micro-PVM" & DF.classifications_0.25_0.09_11676 == "Singlet" & scDblFinder.class == "singlet"), 
#              reduction = "umap.unintegrated", group.by = "DF.classifications_0.25_0.09_11676", label = TRUE, pt.size = 1.5) +
#   labs(title = "Microglia - DoubletFinder Singlets")
p2 = FeaturePlot(merged_obj %>% subset(predicted.subclass == "Micro-PVM"), reduction = "umap.unintegrated", features = "nCount_RNA", order = T, pt.size = 1.5)
print(p1 + p2)
dev.off()

png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/Microglia_markers.png", width = 25, height = 20, units = "in", res = 300)
FeaturePlot(merged_obj %>% subset(predicted.subclass == "Micro-PVM"), 
            features = c("Tmem119", "Aif1", "Cx3cr1", "P2ry12", "Trem2", "Cd3e", "Cd3d", "Cd4", "Cd8a", "Ms4a7", "Cd163", "Mbp"),
            reduction = "umap.unintegrated", order = TRUE, ncol = 4, raster = FALSE)
dev.off()

png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/MS4A_KO_mice_umap_Exclusive.png", width = 8, height = 6, units = "in", res = 300)
DimPlot(merged_obj %>% subset(DF.classifications_0.25_0.09_11676 == "Singlet" & scDblFinder.class == "singlet" ), reduction = "umap.unintegrated", group.by = c("predicted.subclass", ""), label = TRUE) +
  labs(title = "Exclusive") 
dev.off()
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/1-Doublet_removal/MS4A_KO_mice_umap_Inclusive.png", width = 8, height = 6, units = "in", res = 300)
DimPlot(merged_obj %>% subset(DF.classifications_0.25_0.09_11676 == "Singlet" | scDblFinder.class == "singlet" ), reduction = "umap.unintegrated", group.by = c("predicted.subclass", ""), label = TRUE) +
  labs(title = "Inclusive") 
dev.off()

qsave(merged_obj, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/03_MS4A_KO_mice_DoubletFinder_scDblFinder.qs")
merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/03_MS4A_KO_mice_DoubletFinder_scDblFinder.qs")

# merged_obj %>% subset(seurat_clusters %in% c(3, 22)) %>% FindMarkers(ident.1 = "22", ident.2 = "3") -> Cluster3vs22_markers
# merged_obj %>% FindMarkers(ident.1 = "22") -> Cluster22_markers
# Cluster3vs22_markers %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames %>% head(100)
# Cluster22_markers %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames %>% head(100)
# Cluster 22 express microglia markers but also Mbp and inbetween Oligo+Microglia -> might be microglia+Oligo cluster

#########################################################
# DecontX on singlets object
#########################################################

library(SingleCellExperiment) 
library(celda) 

merged_obj = merged_obj %>% subset(DF.classifications_0.25_0.09_11676 == "Singlet" & scDblFinder.class == "singlet")
sce <- as.SingleCellExperiment(merged_obj)
sce <- decontX(sce)
decont_counts <- decontXcounts(sce)

contam <- colData(sce)$decontX_contamination
merged_obj$decontX_contamination <- contam
merged_obj[["decontX"]] <- CreateAssayObject(counts = decont_counts)
DefaultAssay(merged_obj) <- "decontX"

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/2-DecontX")
merged_obj_rm_decont = merged_obj %>% subset(decontX_contamination < 0.2)
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/2-DecontX/MS4A_KO_mice_umap_split_by_celltype.png", width = 20, height = 15, units = "in", res = 300)
DimPlot(merged_obj_rm_decont, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters", "predicted.class","predicted.subclass"), label = TRUE, pt.size = 1.5) + labs(title = "Split by Major Cell Types")
dev.off()

merged_obj_rm_decont = merged_obj
merged_obj_rm_decont$decontX_remove = ifelse(merged_obj_rm_decont$decontX_contamination < 0.2, "retain", "remove")
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/2-DecontX/MS4A_KO_mice_umap_decontx_remove_0.2.png", width = 12, height = 10, units = "in", res = 300)
DimPlot(merged_obj_rm_decont, reduction = "umap.unintegrated", group.by = c("decontX_remove"), label = TRUE, pt.size = 1.5) + labs(title = "Split by Major Cell Types")
dev.off()

merged_obj_rm_decont$decontX_remove = ifelse(merged_obj_rm_decont$decontX_contamination < 0.3, "retain", "remove")
png("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/2-DecontX/MS4A_KO_mice_umap_decontx_remove_0.3.png", width = 12, height = 10, units = "in", res = 300)
DimPlot(merged_obj_rm_decont, reduction = "umap.unintegrated", group.by = c("decontX_remove"), label = TRUE, pt.size = 1.5) + labs(title = "Split by Major Cell Types")
dev.off()

merged_obj = merged_obj %>% subset(decontX_contamination < 0.2)

qsave(merged_obj, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")

# Create a new metadata column for manual annotation
merged_obj$manual_annotation <- "Unknown"  # Initialize with default value

# Assign cell type labels based on seurat_clusters
merged_obj$manual_annotation[merged_obj$seurat_clusters == 0] <- "Oligo"
merged_obj$manual_annotation[merged_obj$seurat_clusters == 3] <- "Micro"
merged_obj$manual_annotation[merged_obj$seurat_clusters == 5] <- "Astro"
merged_obj$manual_annotation[merged_obj$seurat_clusters == 10] <- "OPC"
merged_obj$manual_annotation[merged_obj$seurat_clusters == 14] <- "Other"

# For excitatory neurons (multiple clusters + predicted class filter)
excit_mask <- merged_obj$seurat_clusters %in% c(1, 16, 4, 9, 7, 11, 6, 2, 15, 21, 19, 13, 17, 12) & 
              merged_obj$predicted.class == "Glutamatergic"
merged_obj$manual_annotation[excit_mask] <- "Excit"

# For inhibitory neurons (multiple clusters + predicted class filter)
inhib_mask <- merged_obj$seurat_clusters %in% c(8, 22, 13, 12, 20) & 
              merged_obj$predicted.class == "GABAergic"
merged_obj$manual_annotation[inhib_mask] <- "Inhib"

qsave(merged_obj, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")
merged_obj = qread("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/0-Whole_brain/04_MS4A_KO_mice_decontX.qs")


#########################################################
# Split the object by major cell types
#########################################################
merged_obj@meta.data %>% head

setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/3-Split_by_cell_type")
png("MS4A_KO_mice_umap_split_by_celltype.png", width = 20, height = 15, units = "in", res = 300)
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters", "predicted.class","predicted.subclass"), label = TRUE, pt.size = 1.5) + labs(title = "Split by Major Cell Types")
dev.off()

png("MS4A_KO_mice_umap_split_by_treatment.png", width = 20, height = 15, units = "in", res = 300)
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("orig.ident","predicted.subclass", "genotype", "sex"), label = TRUE, pt.size = 1.5) + labs(title = "Split by Major Cell Types")
dev.off()

Oligo = merged_obj %>% subset(seurat_clusters == 0)
Micro = merged_obj %>% subset(seurat_clusters == 3)
Astro = merged_obj %>% subset(seurat_clusters == 5)
OPC = merged_obj %>% subset(seurat_clusters == 10)
Excit = merged_obj %>% subset(seurat_clusters %in% c(1, 16, 4, 9, 7, 11, 6, 2, 15, 21, 19, 13, 17, 12)) %>% subset(predicted.class == "Glutamatergic") 
Excit@meta.data$predicted.class %>% table
Excit@meta.data$seurat_clusters %>% table
Excit@meta.data$predicted.subclass %>% table
Inhib = merged_obj %>% subset(seurat_clusters %in% c(8, 22, 13, 12, 20)) %>% subset(predicted.class == "GABAergic")
Inhib@meta.data$predicted.class %>% table
Inhib@meta.data$seurat_clusters %>% table
Inhib@meta.data$predicted.subclass %>% table
Other = merged_obj %>% subset(seurat_clusters == 14)


# Harmony Integration for each of the cell types
library(harmony)

#Micro$batch = ifelse(str_detect(Micro$orig.ident, "Exp"), "batch2", "batch1")

Celltype_HarmonyIntegration <- function(obj){
  obj[["decontX"]] <- split(obj[["decontX"]], f = obj$orig.ident)
  obj = obj %>% NormalizeData() %>% FindVariableFeatures()

  # Not use ribosomal gene for variable feature selection
  variable_genes <- VariableFeatures(obj)
  filtered_genes <- variable_genes[!grepl("^Rpl|^Rps", variable_genes)]
  VariableFeatures(obj) <- filtered_genes
  print(table(grepl("^Rpl|^Rps", filtered_genes)))  # Should return all FALSE

  #obj = obj %>% ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) %>% RunPCA()
  obj = obj %>% ScaleData() %>% RunPCA()
  #obj = obj %>% FindNeighbors() 
  #obj = obj %>% FindClusters(resolution = 0.2) 
  #obj = obj %>% RunUMAP(dims = 1:30)
  
  setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/3-Split_by_cell_type")
  png(paste0("MS4A_KO_mice_Elbow_plot_", obj$predicted.subclass[1], ".png"), width = 5, height = 5, units = "in", res = 300)
  p = ElbowPlot(obj, ndims = 50) + labs(title = paste("Elbow Plot for", obj$predicted.subclass[1]))
  print(p)
  dev.off()
  obj = obj %>% IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") 
  obj = JoinLayers(obj)
  return(obj)
}

Celltype_HarmonyIntegration_for_neuron <- function(obj){
  obj[["decontX"]] <- split(obj[["decontX"]], f = obj$orig.ident)
  obj = obj %>% NormalizeData() %>% FindVariableFeatures()

  # Not use ribosomal gene for variable feature selection
  variable_genes <- VariableFeatures(obj)
  filtered_genes <- variable_genes[!grepl("^Rpl|^Rps", variable_genes)]
  VariableFeatures(obj) <- filtered_genes
  print(table(grepl("^Rpl|^Rps", filtered_genes)))  # Should return all FALSE

  #obj = obj %>% ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA")) %>% RunPCA()
  obj = obj %>% ScaleData() %>% RunPCA()
  #obj = obj %>% FindNeighbors() 
  #obj = obj %>% FindClusters(resolution = 0.2) 
  #obj = obj %>% RunUMAP(dims = 1:30)
  
  setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/3-Split_by_cell_type")
  png(paste0("MS4A_KO_mice_Elbow_plot_", obj$predicted.class[1], ".png"), width = 5, height = 5, units = "in", res = 300)
  p = ElbowPlot(obj, ndims = 50) + labs(title = paste("Elbow Plot for", obj$predicted.subclass[1]))
  print(p)
  dev.off()
  obj = obj %>% IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") 
  obj = JoinLayers(obj)
  return(obj)
}


try_ndims = function(obj, range){
  for(ndims in range) {
    obj = obj %>% FindNeighbors(dims = 1:ndims, reduction = "harmony") %>% FindClusters(resolution = 0.2, cluster.name = "harmony_clusters") %>% 
    RunUMAP(dims = 1:ndims, reduction = "harmony", reduction.name = "umap.harmony")
    setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/3-Split_by_cell_type")
    png(paste0("MS4A_KO_mice_umap_harmony_", obj$predicted.subclass[1], "_ndims_", ndims,".png"), width = 20, height = 15, units = "in", res = 300)
    p = DimPlot(obj, reduction = "umap.harmony", group.by = c("orig.ident", "harmony_clusters",  "genotype", "predicted.subclass"), label = TRUE, pt.size = 1.5) + 
    labs(title = paste0("Harmony Integration ", ndims))
    print(p)
    dev.off()
    # png(paste0("MS4A_KO_mice_umap_harmony_", obj$predicted.subclass[1], "_ndims_", ndims,"_by_donor.png"), width = 20, height = 15, units = "in", res = 300)
    # p = DimPlot(obj, reduction = "umap.harmony", group.by = "orig.ident", split.by = "orig.ident", ncol = 4) 
    # print(p)
    # dev.off()
  }
  return(obj)
}

try_ndims_for_neuron = function(obj, range){
  for(ndims in range) {
    obj = obj %>% FindNeighbors(dims = 1:ndims, reduction = "harmony") %>% FindClusters(resolution = 0.2, cluster.name = "harmony_clusters") %>% 
    RunUMAP(dims = 1:ndims, reduction = "harmony", reduction.name = "umap.harmony")
    setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/3-Split_by_cell_type")
    png(paste0("MS4A_KO_mice_umap_harmony_", obj$predicted.class[1], "_ndims_", ndims,".png"), width = 20, height = 15, units = "in", res = 300)
    p = DimPlot(obj, reduction = "umap.harmony", group.by = c("orig.ident", "harmony_clusters",  "genotype", "predicted.subclass"), label = TRUE, pt.size = 1, raster = T) + 
    labs(title = paste0("Harmony Integration ", ndims))
    print(p)
    dev.off()
    # png(paste0("MS4A_KO_mice_umap_harmony_", obj$predicted.subclass[1], "_ndims_", ndims,"_by_donor.png"), width = 20, height = 15, units = "in", res = 300)
    # p = DimPlot(obj, reduction = "umap.harmony", group.by = "orig.ident", split.by = "orig.ident", ncol = 4) 
    # print(p)
    # dev.off()
  }
  return(obj)
}


#Oligodendrocyte
Oligo_Integrated = Celltype_HarmonyIntegration(Oligo)
try_ndims(Oligo_Integrated, range = c(7, 9, 10, 12, 20, 30, 40))
Oligo_Integrated = try_ndims(Oligo_Integrated, range = 40)
qsave(Oligo_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Oligodendrocyte_HarmonyIntegrated.qs")

# Microglia
Micro_Integrated = Celltype_HarmonyIntegration(Micro)
try_ndims(Micro_Integrated, range = c(6, 9, 10, 13, 20, 30, 40)) 
Micro_Integrated = try_ndims(Micro_Integrated, range = 30)
qsave(Micro_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Microglia_HarmonyIntegrated.qs")

#Astrocytes
Astro_Integrated = Celltype_HarmonyIntegration(Astro)
try_ndims(Astro_Integrated, range = c(9, 11, 15, 20, 30, 40))
Astro_Integrated = try_ndims(Astro_Integrated, range = 40)
qsave(Astro_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Astrocytes_HarmonyIntegrated.qs")

#OPC
OPC_Integrated = Celltype_HarmonyIntegration(OPC)
try_ndims(OPC_Integrated, range = c(7, 10, 20, 30, 40))
OPC_Integrated = try_ndims(OPC_Integrated, range = 7)
qsave(OPC_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_OPC_HarmonyIntegrated.qs")

#Excitatory neurons
Excit_Integrated = Celltype_HarmonyIntegration_for_neuron(Excit) #No regress to either features -> cause error in ScaleData
try_ndims_for_neuron(Excit_Integrated, range = c(8, 16, 22, 30, 40))
Excit_Integrated = try_ndims_for_neuron(Excit_Integrated, range = 8)
qsave(Excit_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Excitatory_neurons_HarmonyIntegrated.qs")

# Inhibitory neurons
Inhib_Integrated = Celltype_HarmonyIntegration_for_neuron(Inhib)
try_ndims_for_neuron(Inhib_Integrated, range = c(8, 10, 12, 20, 30, 40))
Inhib_Integrated = try_ndims_for_neuron(Inhib_Integrated, range = 8)
qsave(Inhib_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Inhibitory_neurons_HarmonyIntegrated.qs")

# Other cell types (e.g. VLMC, Pericytes, etc.)
Other_Integrated = Celltype_HarmonyIntegration(Other)
try_ndims(Other_Integrated, range = c(6, 8, 11, 13, 20, 30, 40))
Other_Integrated = try_ndims(Other_Integrated, range = 13)
qsave(Other_Integrated, "/home/yous/Projects/2025-4-5 Ms4a4a KO mice/0-Data/1-Cell_types/MS4A_KO_mice_Other_HarmonyIntegrated.qs")

****************************************
# The regress ribosome and not use ribosomal genes for variable feature selection
# Helps with the integration but still cannot completely remove the batch effect.
# Even integration the batch effect is not completely removed, especially for the Microglia.
# Therefore remove the Exp12 samples from the dataset.
****************************************


#########################################################
# Split the object by major cell types and remove Exp12 samples
#########################################################

#remove Exp12
# merged_obj$batch = ifelse(str_detect(merged_obj$orig.ident, "Exp"), "Exp2", "Exp1")
# merged_obj@meta.data[,c("orig.ident", "batch")] %>% unique
# merged_obj = subset(merged_obj, batch == "Exp1") # keep only Exp1 samples

# setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/4-Split_by_cell_type_remove_Exp12")
# png("MS4A_KO_mice_umap_split_by_celltype.png", width = 20, height = 15, units = "in", res = 300)
# DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters", "predicted.class","predicted.subclass"), label = TRUE, pt.size = 1, raster = T) + labs(title = "Split by Major Cell Types")
# dev.off()
# png("MS4A_KO_mice_umap_split_by_treatment.png", width = 20, height = 15, units = "in", res = 300)
# DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("orig.ident","predicted.subclass", "genotype", "sex"), label = TRUE, pt.size = 1, raster = T) + labs(title = "Split by Major Cell Types")
# dev.off()

# Oligo = merged_obj %>% subset(seurat_clusters == 0)
# Micro = merged_obj %>% subset(seurat_clusters == 3)
# Astro = merged_obj %>% subset(seurat_clusters == 5)
# OPC = merged_obj %>% subset(seurat_clusters == 10)
# Excit = merged_obj %>% subset(predicted.class == "Glutamatergic") %>% subset(seurat_clusters %in% c(1, 16, 4, 9, 7, 11, 6, 2, 15, 21, 19, 13, 17, 12))
# Inhib = merged_obj %>% subset(predicted.class == "GABAergic") %>% subset(seurat_clusters %in% c(8, 22, 13, 12, 20))
# Other = merged_obj %>% subset(seurat_clusters == 14)

# # Harmony Integration for each of the cell types
# library(harmony)

# Celltype_HarmonyIntegration <- function(obj, n){
#   obj[["decontX"]] <- split(obj[["decontX"]], f = obj$orig.ident)
#   obj = obj %>% NormalizeData() %>% FindVariableFeatures()

#   # Not use ribosomal gene for variable feature selection
#   variable_genes <- VariableFeatures(obj)
#   filtered_genes <- variable_genes[!grepl("^Rpl|^Rps", variable_genes)]
#   VariableFeatures(obj) <- filtered_genes
#   print(table(grepl("^Rpl|^Rps", filtered_genes)))  # Should return all FALSE

#   obj = obj %>% ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) %>% RunPCA() 
#   #obj = obj %>% FindNeighbors() 
#   #obj = obj %>% FindClusters(resolution = 0.2) 
#   #obj = obj %>% RunUMAP(dims = 1:30)
  
#   setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/4-Split_by_cell_type_remove_Exp12/")
#   png(paste0("MS4A_KO_mice_Elbow_plot_", obj$predicted.subclass[1], ".png"), width = 5, height = 5, units = "in", res = 300)
#   p = ElbowPlot(obj, ndims = 50) + labs(title = paste("Elbow Plot for", obj$predicted.subclass[1]))
#   print(p)
#   dev.off()
#   obj = obj %>% IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
#   obj = JoinLayers(obj)
#   return(obj)
# }

# try_ndims = function(obj, range){
#   for(ndims in range) {
#     obj = obj %>% FindNeighbors(dims = 1:ndims, reduction = "harmony") %>% FindClusters(resolution = 0.2, cluster.name = "harmony_clusters") %>% 
#     RunUMAP(dims = 1:ndims, reduction = "harmony", reduction.name = "umap.harmony")
#     setwd("/home/yous/Projects/2025-4-5 Ms4a4a KO mice/1-Result/0-Preprocessing/4-Split_by_cell_type_remove_Exp12/")
#     png(paste0("MS4A_KO_mice_umap_harmony_", obj$predicted.subclass[1], "_ndims_", ndims,".png"), width = 20, height = 15, units = "in", res = 300)
#     p = DimPlot(obj, reduction = "umap.harmony", group.by = c("orig.ident", "harmony_clusters",  "genotype", "sex"), label = TRUE, pt.size = 1.5) 
#     print(p)
#     dev.off()
#   }
#   return(obj)
# }

# Oligo_Integrated = Celltype_HarmonyIntegration(Oligo)
# Micro_Integrated = Celltype_HarmonyIntegration(Micro)
# Astro_Integrated = Celltype_HarmonyIntegration(Astro)
# OPC_Integrated = Celltype_HarmonyIntegration(OPC)
# Excit_Integrated = Celltype_HarmonyIntegration(Excit)
# Inhib_Integrated = Celltype_HarmonyIntegration(Inhib)
# Other_Integrated = Celltype_HarmonyIntegration(Other)

# try_ndims(Oligo_Integrated, range = c(7, 12, 20, 30, 40))
# try_ndims(Micro_Integrated, range = c(7, 12, 20, 30, 40))

# try_ndims(Astro_Integrated, range = c(6, 9, 10, 13, 20, 30))
# try_ndims(OPC_Integrated, range = c(6, 9, 10, 13, 20, 30))
# try_ndims(Excit_Integrated, range = c(6, 9, 10, 13, 20, 30))
# try_ndims(Inhib_Integrated, range = c(6, 9, 10, 13, 20, 30))
# try_ndims(Other_Integrated, range = c(6, 9, 10, 13, 20, 30)) 

# ****************************************
# Later found Exp12 after remove batch effect, they are acting like more homeostatic microglia, so retain these cells for analysis.
# ****************************************


