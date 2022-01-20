library(dplyr)
library(Seurat)

#Create Vector of Cell Type Labels and Extract Gene Expression Matrix from Seurat Object
Adipose = Read10X("SVF")
Adipose2 = CreateSeuratObject(Adipose)
Adipose3 <- subset(Adipose2, downsample = 1000)
write.table(Adipose3@assays[["RNA"]]@counts, file='Gene_Counts_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)
metadata = read.delim("SVF/metaData.tsv")
metadata = metadata[,-1]
write.csv(metadata,"Adipose_Cell_Labels.csv")
write.table(Adipose,"Gene_Counts_Per_Cell.tsv")

#Combine cell type labels and Gene Expression Matrix in Python



#Add Genotype data for each sample to Cibersortx Results
data = read.csv("CIBERSORTx_Job31_Results.csv")
genotypes = read.csv("Adipose_3wa_raw&processed_2/expr.lumi.gene(55_9400)_vst_0.01_dr0.3_3wayANOVA_T.csv")
 genotypes = genotypes[9401:9404,]
 genotypes = t(genotypes)
 colnames(genotypes)= genotypes[1,]
 genotypes = genotypes[-1,]
 genotypes = genotypes[-(37:50),]
 final = cbind(genotypes,data)
 write.csv(final, "Adipose_Cibersort_3w_T2.csv")
