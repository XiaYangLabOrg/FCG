#=========ING==========
ing <- read.table("ING/expr.lumi.gene(55_9400)_vst_0.01_dr0.3.txt", stringsAsFactors = F)
ing <- t(ing)
ing <- as.data.frame(ing)
sample.info <- as.data.frame(read.table("ING/Group_info_ING.txt", header = T))
ing$group <- sample.info$Sample.Group[match(rownames(ing), sample.info$Sample.ID)]

ing$chr[ing$group %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
ing$chr[ing$group %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY-"
ing$chr <- factor(ing$chr)

ing$sry[ing$group %in% c("XY-Sry_B", "XY-Sry_E", "XY-Sry_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "Sry"
ing$sry[ing$group %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "-"
ing$sry <- factor(ing$sry)

ing$treatment[ing$group %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
ing$treatment[ing$group %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
ing$treatment[ing$group %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
ing$treatment <- factor(ing$treatment)

ing[1:5,9400:9404]
chr_g <- unique(ing$chr)
sry_g <- unique(ing$sry)
treatment_g <- unique(ing$treatment)

group <- as.data.frame(ing[1:55,9401:9404])
library(plyr)
group[35:40,]
group$chr <- revalue(group$chr, c("XY-" = "XY"))
group$sry <- revalue(group$sry, c("-" = "NULL"))
chr_g <- unique(group$chr)
sry_g <- unique(group$sry)
treatment_g <- unique(group$treatment)

ing <- t(ing[1:55, 1:9400])
ing <- as.data.frame(ing)
# genes <- rownames(ing)
# ing <- as.data.frame(apply(ing, 2, 
#                            function(x){as.numeric(as.character(x))}))
# rownames(ing) <- genes
ing[1:5,1:5]
library(limma)
for(i in 1:2){
  for(j in 1:2){
    chr <- as.character(chr_g[i])
    sry <- as.character(sry_g[j])
    ing_sub <- ing[,which(group$chr == chr & group$sry == sry)]
    group_sub <- group[which(group$chr == chr & group$sry == sry),]
    design = model.matrix(~0 + treatment, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(ing_sub, design = design)
    contrast <- makeContrasts(treatmentT - treatmentE, treatmentE-treatmentB, treatmentT-treatmentB, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("ING/genelist/T-E_0.3_",chr,"_",sry,".xlsx"))
    
    tt2 <- topTable(raw.fit, coef = 2, n = Inf, adjust.method = "BH")
    tt2 <- tt2[order(tt2$adj.P.Val, decreasing = F),]
    tt2$gene <- rownames(tt2)
    openxlsx::write.xlsx(tt2, 
                         file = paste0("ING/genelist/E-B_0.3_",chr,"_",sry,".xlsx"))
    
    tt3 <- topTable(raw.fit, coef = 3, n = Inf, adjust.method = "BH")
    tt3 <- tt3[order(tt3$adj.P.Val, decreasing = F),]
    tt3$gene <- rownames(tt3)
    openxlsx::write.xlsx(tt3, 
                         file = paste0("ING/genelist/T-B_0.3_",chr,"_",sry,".xlsx"))
  }
}
  


for(i in 1:2){
  for(j in 1:3){
    chr <- as.character(chr_g[i])
    treatment <- as.character(treatment_g[j])
    ing_sub <- ing[,which(group$chr == chr & group$treatment == treatment)]
    group_sub <- group[which(group$chr == chr & group$treatment == treatment),]
    design = model.matrix(~0 + sry, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(ing_sub, design = design)
    contrast <- makeContrasts(srySry - sryNULL, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("ING/genelist/Sry-Null_0.3_",chr,"_",treatment,".xlsx"))
  }
}

for(i in 1:2){
  for(j in 1:3){
    sry <- as.character(sry_g[i])
    treatment <- as.character(treatment_g[j])
    ing_sub <- ing[,which(group$sry == sry & group$treatment == treatment)]
    group_sub <- group[which(group$sry == sry & group$treatment == treatment),]
    design = model.matrix(~0 + chr, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(ing_sub, design = design)
    contrast <- makeContrasts(chrXX - chrXY, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("ING/genelist/XX-XY_0.3_new",sry,"_",treatment,".xlsx"))
  }
}

#===========liver===========
liver <- read.table("liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.3.txt", stringsAsFactors = F)
liver <- t(liver)
liver <- as.data.frame(liver)
sample.info <- as.data.frame(read.table("liver/Group_info_liver.txt", header = T))
sample.info$Sample.ID <- paste0("X", sample.info$Sample.ID)
liver$group <- sample.info$Sample.Group[match(rownames(liver), sample.info$Sample.ID)]

liver$chr[liver$group %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
liver$chr[liver$group %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY-"
liver$chr <- factor(liver$chr)

liver$sry[liver$group %in% c("XY-Sry_B", "XY-Sry_E", "XY-Sry_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "Sry"
liver$sry[liver$group %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "-"
liver$sry <- factor(liver$sry)

liver$treatment[liver$group %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
liver$treatment[liver$group %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
liver$treatment[liver$group %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
liver$treatment <- factor(liver$treatment)

liver[1:5,7466:7470]

group <- as.data.frame(liver[1:60,7467:7470])
library(plyr)
group[35:40,]
group$chr <- revalue(group$chr, c("XY-" = "XY"))
group$sry <- revalue(group$sry, c("-" = "NULL"))
chr_g <- unique(group$chr)
sry_g <- unique(group$sry)
treatment_g <- unique(group$treatment)

liver <- t(liver[1:60, 1:7466])
liver <- as.data.frame(liver)
# genes <- rownames(liver)
# liver <- as.data.frame(apply(liver, 2, 
#                            function(x){as.numeric(as.character(x))}))
# rownames(liver) <- genes
liver[1:5,1:5]
library(limma)
for(i in 1:2){
  for(j in 1:2){
    chr <- as.character(chr_g[i])
    sry <- as.character(sry_g[j])
    liver_sub <- liver[,which(group$chr == chr & group$sry == sry)]
    group_sub <- group[which(group$chr == chr & group$sry == sry),]
    design = model.matrix(~0 + treatment, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(liver_sub, design = design)
    contrast <- makeContrasts(treatmentT - treatmentE, treatmentE-treatmentB, treatmentT-treatmentB, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("liver/genelist/T-E_0.3_",chr,"_",sry,".xlsx"))

    tt2 <- topTable(raw.fit, coef = 2, n = Inf, adjust.method = "BH")
    tt2 <- tt2[order(tt2$adj.P.Val, decreasing = F),]
    tt2$gene <- rownames(tt2)
    openxlsx::write.xlsx(tt2, 
                         file = paste0("liver/genelist/E-B_0.3_",chr,"_",sry,".xlsx"))

    tt3 <- topTable(raw.fit, coef = 3, n = Inf, adjust.method = "BH")
    tt3 <- tt3[order(tt3$adj.P.Val, decreasing = F),]
    tt3$gene <- rownames(tt3)
    openxlsx::write.xlsx(tt3, 
                         file = paste0("liver/genelist/T-B_0.3_",chr,"_",sry,".xlsx"))
  }
  
}

for(i in 1:2){
  for(j in 1:3){
    chr <- as.character(chr_g[i])
    treatment <- as.character(treatment_g[j])
    liver_sub <- liver[,which(group$chr == chr & group$treatment == treatment)]
    group_sub <- group[which(group$chr == chr & group$treatment == treatment),]
    design = model.matrix(~0 + sry, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(liver_sub, design = design)
    contrast <- makeContrasts(srySry - sryNULL, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("liver/genelist/Sry-Null_0.3_",chr,"_",treatment,".xlsx"))
  }
}

for(i in 1:2){
  for(j in 1:3){
    sry <- as.character(sry_g[i])
    treatment <- as.character(treatment_g[j])
    liver_sub <- liver[,which(group$sry == sry & group$treatment == treatment)]
    group_sub <- group[which(group$sry == sry & group$treatment == treatment),]
    design = model.matrix(~0 + chr, data = group_sub)
    design[1:5,]
    raw.fit <- limma::lmFit(liver_sub, design = design)
    contrast <- makeContrasts(chrXX - chrXY, levels = design)
    raw.fit <- contrasts.fit(raw.fit, contrasts = contrast)
    raw.fit <- eBayes(raw.fit, trend = T, robust = T)
    tt1 <- topTable(raw.fit, coef = 1, n = Inf, adjust.method = "BH")
    tt1 <- tt1[order(tt1$adj.P.Val, decreasing = F),]
    tt1$gene <- rownames(tt1)
    openxlsx::write.xlsx(tt1, 
                         file = paste0("liver/genelist/XX-XY_0.3_",sry,"_",treatment,".xlsx"))
  }
}
