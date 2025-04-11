# require("biomaRt")
# mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# convert_ensg_to_symbol <- function(ensg_list) {
#     gene_symbols <- getBM(
#         filters = "ensembl_gene_id",
#         attributes = c("ensembl_gene_id", "hgnc_symbol"),
#         values = ensg_list,
#         mart = mart
#     )
#     merged_result <- merge(data.frame(ensembl_gene_id = ensg_list), gene_symbols, by = "ensembl_gene_id", all.x = TRUE)
#     return(merged_result$hgnc_symbol)
# }

# Example usage
# ensg_list <- c("ENSG00000139618", "ENSG00000227232")
# gene_symbols <- convert_ensg_to_symbol(ensg_list)
# print(gene_symbols)


require("org.Hs.eg.db")
require("AnnotationDbi")

convert_ensg_to_symbol_bioconductor <- function(ensg_list) {
    gene_symbols <- mapIds(org.Hs.eg.db, keys = ensg_list, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    return(gene_symbols)
}

# Example usage
# ensg_list <- c("ENSG00000139618", "ENSG00000227232")
# gene_symbols <- convert_ensg_to_symbol_bioconductor(ensg_list)
# print(gene_symbols)