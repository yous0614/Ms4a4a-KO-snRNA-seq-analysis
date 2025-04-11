##########################################################
# Author: Shih-Feng You
# Affiliation: Karch Lab, Washington University in St. Louis, St. Louis, MO, USA
# Date: 2025-4-6
# Description: A function to create a volcano plot
##########################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

date = Sys.Date()

# Define a general-purpose volcano plot function
volcano_plot <- function(data, pval_col, fc_col, path, title, gene = NULL, genelist = NULL, pval_threshold = 0.05, fc_threshold = 0, wid = 8, hei = 6, x = NULL, y = NULL) {
  
  # Ensure the gene column exists in data
  # if (!gene %in% colnames(data)) {
  #   stop("The specified gene column does not exist in the data.")
  # }
  
  # Infer gene names from rownames if gene column not provided
  if (is.null(gene)) {
    data$gene <- rownames(data)
    gene <- "gene"
  }
  
  # Filter out unwanted genes
  data <- data[str_subset(data[[gene]], "Rik$|^Gm|^mt-", negate = TRUE), ]
  
  #Add artificial limit
  if(!is.null(x)){
    data[[fc_col]][data[[fc_col]] > x] <- x
    data[[fc_col]][data[[fc_col]] < -x] <- -x
  }
  if(!is.null(y)){
    data[[pval_col]][-log(data[[pval_col]]) > y] <- -log(y)
  }

  # Create a new column for significance based on thresholds
  data$Significance <- ifelse(data[[pval_col]] < pval_threshold & abs(data[[fc_col]]) > fc_threshold, 
                              ifelse(data[[fc_col]] > fc_threshold, "Upregulated", "Downregulated"),
                              "Not Significant")
  
  # Create a label column for genes in genelist or significant genes if genelist is missing
  data$Label <- if(is.null(genelist)) {
    ifelse(data[[pval_col]] < pval_threshold & abs(data[[fc_col]]) > fc_threshold, data[[gene]], NA)
  } else {
    ifelse(data[[gene]] %in% genelist, data[[gene]], NA)
  }
  
  # Plot using ggplot
  png(paste0(path, date, "_Volcano_", title, ".png"), width = wid, height = hei, res = 300, units = "in")
  p = ggplot(data, aes_string(x = fc_col, y = paste0("-log10(", pval_col, ")"), color = "Significance")) +
    geom_point() +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "royal blue", "Not Significant" = "grey")) +
    geom_text_repel(data = subset(data, !is.na(Label)), aes(label = Label), segment.alpha = 0.5, direction = "both", size = 4.5, box.padding = 0.5, max.overlaps = 20, bg.color = alpha(c("white"), 1), bg.r = 0.3) +
    labs(x = paste0(fc_col), y = paste0("-log10(", pval_col, ")"), title = title) +
    theme_classic2()
  print(p)
  dev.off()
}


#Example usage
# volcano_plot(data = data, pval_col = "Adjusted.P.value", fc_col = "log2FoldChange", gene = "Gene", path = "./", title = "Volcano_Plot", genelist = c("Gene1", "Gene2"), pval_threshold = 0.05, fc_threshold = 0)


volcano_plot_two_threshold <- function(data, pval_col, padj_col, fc_col, gene, path, title, genelist = NULL, pval_threshold = 0.05, fc_threshold = 0, wid = 8, hei = 6, x = NULL, y = NULL) {
  
  # Ensure the gene column exists in data
  if (!gene %in% colnames(data)) {
    stop("The specified gene column does not exist in the data.")
  }
  
  #Add artificial limit
  if(!is.null(x)){
    data[[fc_col]][data[[fc_col]] > x] <- x
    data[[fc_col]][data[[fc_col]] < -x] <- -x
  }
  if(!is.null(y)){
    data[[pval_col]][data[[pval_col]] > y] <- y
  }

  # Create a new column for significance based on thresholds
  data$Significance <- ifelse(data[[pval_col]] < pval_threshold & abs(data[[fc_col]]) > fc_threshold, 
                              ifelse(data[[fc_col]] > fc_threshold, "Upregulated (raw p)", "Downregulated (raw p)"),
                              "Not Significant")

  # Create a new column for significance based on adjusted p-value
  data$Significance <- ifelse(data[[padj_col]] < pval_threshold & abs(data[[fc_col]]) > fc_threshold, 
                              ifelse(data[[fc_col]] > fc_threshold, "Upregulated (padj)", "Downregulated (padj)"),
                              data$Significance)                
  
  # Create a label column for genes in genelist or significant genes if genelist is missing
  data$Label <- if(is.null(genelist)) {
    ifelse(data[[pval_col]] < pval_threshold & abs(data[[fc_col]]) > fc_threshold, data[[gene]], NA)
  } else {
    ifelse(data[[gene]] %in% genelist, data[[gene]], NA)
  }
  
  # Plot using ggplot
  png(paste0(path, date, "_Volcano_", title, ".png"), width = wid, height = hei, res = 300, units = "in")
  p = ggplot(data, aes_string(x = fc_col, y = paste0("-log10(", pval_col, ")"), color = "Significance")) +
    geom_point() +
    scale_color_manual(values = c("Upregulated (padj)" = "red", "Upregulated (raw p)" = "pink", "Downregulated (padj)" = "royal blue", "Downregulated (raw p)" = "#6395EE", "Not Significant" = "grey")) +
    geom_text_repel(data = subset(data, !is.na(Label)), aes(label = Label), segment.alpha = 0.5, direction = "both", size = 4.5, box.padding = 0.5, max.overlaps = 20, bg.color = alpha(c("white"), 1), bg.r = 0.3) +
    labs(x = paste0(fc_col), y = paste0("-log10(", pval_col, ")"), title = title) +
    theme_classic2()
  print(p)
  dev.off()
}


alternative_volcano_plot <- function(data, pval_col = "p_val_adj", fc_col = "avg_log2FC", gene = NULL,
                                     path = "./", title = "Volcano", color = "black",
                                     pval_thresh = 0.05, xlimit = NULL, wid = 10, hei = 6,
                                     manual_add_genes = NULL, n = 30, comparison = "Comparison", 
                                     genotype = "", celltype = "") {
  # Use system date
  date = Sys.Date()
  
  # Infer gene names from rownames if gene column not provided
  if (is.null(gene)) {
    data$gene <- rownames(data)
    gene <- "gene"
  }
  
  # Filter out unwanted genes
  data <- data[str_subset(data[[gene]], "Rik$|^Gm|^mt-", negate = TRUE), ]
  
  # Filter based on p-value threshold
  data <- data %>% filter(.data[[pval_col]] < pval_thresh)

  # Apply x-axis limit (log2FC)
  if (!is.null(xlimit)) {
    data[[fc_col]] <- ifelse(data[[fc_col]] < -xlimit, -xlimit, data[[fc_col]])
    data[[fc_col]] <- ifelse(data[[fc_col]] > xlimit, xlimit, data[[fc_col]])
  }
  
  # Add mean expression level
  data$MeanExpression <- (data$pct.1 + data$pct.2) / 2

  # Select top genes
  top_genes <- data[order(data[[pval_col]]), ][1:n, ]
  
  # Add manually specified genes
  if (!is.null(manual_add_genes)) {
    manual_genes_df <- data %>% filter(.data[[gene]] %in% manual_add_genes)
    top_genes <- bind_rows(manual_genes_df, top_genes) %>% distinct(.data[[gene]], .keep_all = TRUE)
    print(top_genes)
  }

  # Generate plot
  png(paste0(path, date, "_Volcano_", title, ".png"), width = wid, height = hei, units = "in", res = 300)
  p <- ggplot(data, aes_string(x = fc_col, y = "MeanExpression")) +
    geom_point(aes_string(size = paste0("-log10(", pval_col, ")")), color = color, alpha = 0.6) +
    scale_size_continuous(range = c(0.5, 5)) +  
    geom_text_repel(data = top_genes, aes_string(label = gene, x = fc_col, y = "MeanExpression"), size = 6) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),  
          axis.title = element_text(size = 14),
          legend.position = c(0.95, 0.75),
          legend.justification = c(1, 1)) +
    labs(title = paste0("Ms4a4a ", genotype, " ", celltype, " ", comparison), 
         x = fc_col, y = "Mean Expression Level")

  print(p)
  dev.off()
}

# alternative_volcano_plot(data = your_data,
#                          fc_col = "avg_log2FC",
#                          pval_col = "p_val_adj",
#                          gene = "gene",
#                          path = "./",
#                          title = "My_Alt_Volcano",
#                          color = "firebrick",
#                          genotype = "KO",
#                          celltype = "Microglia",
#                          comparison = "PolyIC vs VEH",
#                          manual_add_genes = c("Trem2", "Apoe"))


# alternative_volcano_plot = function(DEGtable, Celltype, Color, Path, Genotype, xlimit = NULL, wid = 10, hei = 6, 
#                             pval_col = "p_val_adj", pval_thresh = 0.05, manual_add_genes = NULL, n = 30, comparison = "PolyIC vs VEH") {
#   # Ensure the chosen p-value column exists
#   if (!pval_col %in% colnames(DEGtable)) {
#     stop("Specified p-value column does not exist in DEGtable.")
#   }
  
#   DEGtable = DEGtable[str_subset(row.names(DEGtable), "Rik$|^Gm|^mt-", negate = TRUE), ]
  
#   # Filter based on the selected p-value column
#   DEGtable = DEGtable %>% filter(.data[[pval_col]] < pval_thresh)
  
#   if (!is.null(xlimit)) {
#     DEGtable$avg_log2FC = ifelse(DEGtable$avg_log2FC < -xlimit, -xlimit, DEGtable$avg_log2FC)
#     DEGtable$avg_log2FC = ifelse(DEGtable$avg_log2FC > xlimit, xlimit, DEGtable$avg_log2FC)
#   }
  
#   DEGtable$MeanExpression = (DEGtable$pct.1 + DEGtable$pct.2) / 2
#   DEGtable$gene = rownames(DEGtable)
  
#   # Select top genes based on p-value
#   top_genes = DEGtable[order(DEGtable[[pval_col]]), ][1:n, ]
  
#   # Ensure manual_add_genes exist in the dataset
#   if (!is.null(manual_add_genes)) {
#     manual_genes_df = DEGtable %>% filter(gene %in% manual_add_genes)
#     top_genes = bind_rows(manual_genes_df, top_genes) %>% distinct(gene, .keep_all = TRUE)
#     print(top_genes)
#   }

#   png(paste0(Path, "Trem2", Genotype, "_", comparison, "_", Celltype, "_volcano.png"), width = wid, height = hei, units = "in", res = 300)
  
#   p = ggplot(DEGtable, aes(x = avg_log2FC, y = MeanExpression)) +
#     geom_point(aes(size = -log10(.data[[pval_col]])), color = Color, alpha = 0.6) +
#     scale_size_continuous(range = c(0.5, 5)) +  
#     geom_text_repel(data = top_genes, aes(label = gene, x = avg_log2FC, y = MeanExpression), size = 6) +
#     theme_classic() +
#     theme(axis.text = element_text(size = 12),  
#           axis.title = element_text(size = 14),
#           legend.position = c(0.95, 0.75),
#           legend.justification = c(1, 1)) +
#     labs(title = paste0("Trem2", Genotype, " ", Celltype, " ", comparison), x = "log2FC", y = "Mean Expression Level")
  
#   print(p)
#   dev.off()
# }


