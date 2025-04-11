##########################################################
# Author: Shih-Feng You
# Affiliation: Karch Lab, Washington University in St. Louis, St. Louis, MO, USA
# Date: 2025-4-6
# Description: A function to perform EnrichR analysis for DEG tables
##########################################################

library(enrichR)
library(openxlsx)
library(stringr)
library(ggplot2)
library(ggpubr)

# Define a general-purpose EnrichR function
today = Sys.Date()

enrichR_pathway = function(data, pval_col, fc_col, gene, path, title, today = Sys.Date(), databases = c("KEGG_2019_Mouse", "GO_Biological_Process_2023", "GO_Cellular_Component_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "WikiPathways_2024_Mouse")) {
    data = data[!grepl("Rik$|^Gm|^mt-|^Rpl|^Rps", data[[gene]]), ]
    data = data[data[[pval_col]] < 0.05, ]
    upDEG = data[data[[fc_col]] > 0 & data[[pval_col]] < 0.05, gene]
    downDEG = data[data[[fc_col]] < 0 & data[[pval_col]] < 0.05, gene]
    print(paste0("Upregulated DEGs: ", length(upDEG)))
    print(paste0("Downregulated DEGs: ", length(downDEG)))

    sheet = list()
    for(database in databases){
        
        if (length(upDEG) > 0) {
            UpPath = enrichr(upDEG, database)
            UpPath = UpPath[[database]][, c("Term", "P.value", "Adjusted.P.value", "Genes")] %>% filter(P.value < 0.05)
        } else {
            UpPath = data.frame(Term = character(0), P.value = numeric(0), Adjusted.P.value = numeric(0), Genes = character(0))
        }
        
        if (length(downDEG) > 0) {
            DownPath = enrichr(downDEG, database)
            DownPath = DownPath[[database]][, c("Term", "P.value", "Adjusted.P.value", "Genes")] %>% filter(P.value < 0.05)
        } else {
            DownPath = data.frame(Term = character(0), P.value = numeric(0), Adjusted.P.value = numeric(0), Genes = character(0))
        }

        if (nrow(UpPath) > 0) {
            UpPath$Direction = "Up"
        } else {
            UpPath$Direction = character(0)
        }
        if (nrow(DownPath) > 0) {
            DownPath$Direction = "Down"
        } else {
            DownPath$Direction = character(0)
        }
        
        Path = rbind(UpPath, DownPath)
        Path = Path[,c("Direction", "Term", "P.value", "Adjusted.P.value", "Genes")]
        sheet[[database]] = Path
    }
    write.xlsx(sheet, paste0(path, today, "_", title, "_EnrichR.xlsx"), asTable = TRUE)
    return(sheet)
}

enrichR_bar_plot = function(df, title = NULL){
    df = df %>% filter(Adjusted.P.value < 0.05)
    df <- df %>%
        mutate(log_adj_pval = -log10(Adjusted.P.value),
               log_adj_pval = ifelse(Direction == "Up", log_adj_pval, -log_adj_pval))
    df$bar_color <- ifelse(df$Adjusted.P.value > 0.05, "gray80", ifelse(df$Direction == "Up", "red", "blue"))
    p = ggplot(df, aes(x = log_adj_pval, y = reorder(Term, log_adj_pval), fill = bar_color)) +
        geom_bar(stat = "identity") +
        scale_fill_identity() +
        labs(x = "-log10(Adjusted P-value)", y = "Term", title = title) +
        theme_classic2() +
        theme(text = element_text(size = 16), axis.text.y = element_text(size = 20))
    return(p)
}
