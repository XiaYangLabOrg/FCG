### convert mouse genes into human genes in mouse modules

genes <- read.csv("modules_for_GO_ing.csv",header = T)
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  # humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}
genes_H <- convertMouseGeneList(genes$Symbol)
genes$human <- genes_H$HGNC.symbol[match(genes$Symbol, genes_H$MGI.symbol)]
head(genes)
genes$human[which(is.na(genes$human))] <- toupper(genes$Symbol[which(is.na(genes$human))])
write.csv(genes, "modules_for_GO_mapped.csv")