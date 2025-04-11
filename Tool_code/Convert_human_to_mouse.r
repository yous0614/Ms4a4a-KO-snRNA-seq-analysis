
library(dplyr)

# mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
# mouse_human_genes %>% write.csv("/hpfs/userws/yous09/Tool code/mouse_human_genes.csv",row.names = FALSE)
mouse_human_genes = read.csv("/hpfs/userws/yous09/0-Function/Convert_human_and_mouse/mouse_human_genes.csv")

 convert_mouse_to_human <- function(gene_list) { 
      output = c()
      #mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

      for(gene in gene_list) {
           class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
           if( !identical(class_key, integer(0)) ) {
                human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
                for(human_gene in human_genes) {
                     output = rbind(c(gene, human_gene), output)
                }
           }
      }
      return (output[,2] %>% rev)
 }

 convert_human_to_mouse <- function(gene_list) {
     output = c()
     #mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

     for(gene in gene_list) {
           class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
           if( !identical(class_key, integer(0)) ) {
             human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
             for(human_gene in human_genes) {
                 output = rbind(c(gene, human_gene), output)
             }
           }
      }
      return (output[,2] %>% rev)
 }

