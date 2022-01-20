#======liver 1st PC=======
expr_liver <- read.table("liver/liver_expr.lumi.gene(60_7466)_vst_0.01_dr0.3.txt")
expr_liver[1:5,1:5]

moduleTable_liver <- read.csv("moduleTable.csv")
head(moduleTable_liver)
module_id <- as.character(moduleTable_liver$module.id[which(moduleTable_liver$module.size > 50)])

module_liver <- read.csv("modules.csv")
module_liver[1:5,1:5]
rownames(module_liver) <- module_liver$X
module_liver <- module_liver[,-1]
module_liver <- module_liver[module_id,]
module_liver_n <- apply(module_liver, 2, function(x){as.character(x)})
rownames(module_liver_n) <- rownames(module_liver)

network_liver <- matrix(NA, nrow = 92, ncol = 60)
network_liver <- as.data.frame(network_liver)
rownames(network_liver) <- module_id
colnames(network_liver) <- colnames(expr_liver)

for(i in 1:length(module_id)){
  expr_sub <- expr_liver[which(rownames(expr_liver) %in% module_liver_n[i,]),]
  pca <- prcomp(t(expr_sub), center = T)
  network_liver[i,] <- pca$x[,1]
}
write.csv(network_liver, "network_liver.csv")

#=========ING 1st PC==========
expr_ing <- read.table("ING/ING_expr.lumi.gene(55_9400)_vst_0.01_dr0.3.txt")
expr_ing[1:5,1:5]

moduleTable_ing <- read.csv("moduleTable_ing.csv")
head(moduleTable_ing)
module_id <- as.character(moduleTable_ing$module.id[which(moduleTable_ing$module.size > 50)])

module_ing <- read.csv("modules_ing.csv")
module_ing[1:5,1:5]
rownames(module_ing) <- module_ing$X
module_ing <- module_ing[,-1]
module_ing <- module_ing[module_id,]
module_ing_n <- apply(module_ing, 2, function(x){as.character(x)})
rownames(module_ing_n) <- rownames(module_ing)

network_ing <- matrix(NA, nrow = 128, ncol = 55)
network_ing <- as.data.frame(network_ing)
rownames(network_ing) <- module_id
colnames(network_ing) <- colnames(expr_ing)

for(i in 1:length(module_id)){
  expr_sub <- expr_ing[which(rownames(expr_ing) %in% module_ing_n[i,]),]
  pca <- prcomp(t(expr_sub), center = T)
  network_ing[i,] <- pca$x[,1]
}
write.csv(network_ing, "network_ing.csv")

#=====3-way======
sample.info_ing <- as.data.frame(read.table("ING/Group_info_ING.txt", header = T))
ing <- t(network_ing)
ing <- as.data.frame(ing)

ing$group <- sample.info_ing$Sample.Group[match(rownames(ing), sample.info_ing$Sample.ID)]

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

chr.p <- NA
sry.p <- NA
treatment.p <- NA
`chr:sry` <- NA
`chr:treatment` <- NA
`sry:treatment` <- NA
`chr:sry:treatment` <- NA

ing_T <- ing[which(ing$treatment != "E"),]
ing_E <- ing[which(ing$treatment != "T"),]
ing_E <- rbind(ing_E, chr.p, sry.p, treatment.p, `chr:sry`, `chr:treatment`, `sry:treatment`, `chr:sry:treatment`)
ing_T <- rbind(ing_T, chr.p, sry.p, treatment.p, `chr:sry`, `chr:treatment`, `sry:treatment`, `chr:sry:treatment`)
rownames(ing_E)[37:43] <- c("chr.p", "sry.p", "treatment.p", "`chr:sry`", "`chr:treatment`", "`sry:treatment`", "`chr:sry:treatment`")
rownames(ing_T)[37:43] <- c("chr.p", "sry.p", "treatment.p", "`chr:sry`", "`chr:treatment`", "`sry:treatment`", "`chr:sry:treatment`")

for(i in 1:128){
  fm <- aov(unlist(ing_T[1:36,i]) ~ chr * sry * treatment, data = ing_T[1:36,])
  ing_T[37:43, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}

for(i in 1:128){
  fm <- aov(unlist(ing_E[1:36,i]) ~ chr * sry * treatment, data = ing_E[1:36,])
  ing_E[37:43, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}

liver <- t(network_liver)
liver <- as.data.frame(liver)
sample.info <- as.data.frame(read.table("liver/Group_info_liver.txt", header = T))
sample.info$Sample.ID <- paste0("X", sample.info$Sample.ID)
# sample.info <- sample.info[-which(sample.info$Sample.ID %in% bad),]
# 
# write.table(sample.info, "liver/Targets.new.txt", row.names = F, col.names = T)

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

chr.p <- NA
sry.p <- NA
treatment.p <- NA
`chr:sry` <- NA
`chr:treatment` <- NA
`sry:treatment` <- NA
`chr:sry:treatment` <- NA

liver_E <- liver[which(liver$treatment != "T"),]
liver_T <- liver[which(liver$treatment != "E"),]
liver_E <- rbind(liver_E, chr.p, sry.p, treatment.p, `chr:sry`, `chr:treatment`, `sry:treatment`, `chr:sry:treatment`)
liver_T <- rbind(liver_T, chr.p, sry.p, treatment.p, `chr:sry`, `chr:treatment`, `sry:treatment`, `chr:sry:treatment`)
rownames(liver_T)[41:47] <- c("chr.p", "sry.p", "treatment.p", "`chr:sry`", "`chr:treatment`", "`sry:treatment`", "`chr:sry:treatment`")
rownames(liver_E)[41:47] <- c("chr.p", "sry.p", "treatment.p", "`chr:sry`", "`chr:treatment`", "`sry:treatment`", "`chr:sry:treatment`")

for(i in 1:92){
  fm <- aov(unlist(liver_T[1:40,i]) ~ chr * sry * treatment, data = liver_T[1:40,])
  liver_T[41:47, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}
for(i in 1:92){
  fm <- aov(unlist(liver_E[1:40,i]) ~ chr * sry * treatment, data = liver_E[1:40,])
  liver_E[41:47, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}


ing_T <- as.data.frame(t(ing_T))
ing_E <- as.data.frame(t(ing_E))
liver_T <- as.data.frame(t(liver_T))
liver_E <- as.data.frame(t(liver_E))

ing_E$chr.p.adj <- NA
ing_E$sry.p.adj <- NA
ing_E$treatment.p.adj <- NA
ing_E$chr.sry.p.adj <- NA
ing_E$chr.treatment.p.adj <- NA
ing_E$sry.treatment.p.adj <- NA
ing_E$chr.sry.treatment.p.adj <- NA

ing_T$chr.p.adj <- NA
ing_T$sry.p.adj <- NA
ing_T$treatment.p.adj <- NA
ing_T$chr.sry.p.adj <- NA
ing_T$chr.treatment.p.adj <- NA
ing_T$sry.treatment.p.adj <- NA
ing_T$chr.sry.treatment.p.adj <- NA

ing_T[,37:43] <- apply(ing_T[,37:43], 2, function(x){as.numeric(as.character(x))})
ing_E[,37:43] <- apply(ing_E[,37:43], 2, function(x){as.numeric(as.character(x))})

for(i in 37:43){
  ing_E[1:128,i+7] <- p.adjust(ing_E[1:128,i], method = 'fdr')
}
for(i in 37:43){
  ing_T[1:128,i+7] <- p.adjust(ing_T[1:128,i], method = 'fdr')
}
ing_E[128:132, 37:50]

write.csv(ing_T, "ING/network_ing_3way_T.csv", row.names = T, col.names = T)
write.csv(ing_E, "ING/network_ing_3way_E.csv", row.names = T, col.names = T)

liver_E$chr.p.adj <- NA
liver_E$sry.p.adj <- NA
liver_E$treatment.p.adj <- NA
liver_E$chr.sry.p.adj <- NA
liver_E$chr.treatment.p.adj <- NA
liver_E$sry.treatment.p.adj <- NA
liver_E$chr.sry.treatment.p.adj <- NA

liver_T$chr.p.adj <- NA
liver_T$sry.p.adj <- NA
liver_T$treatment.p.adj <- NA
liver_T$chr.sry.p.adj <- NA
liver_T$chr.treatment.p.adj <- NA
liver_T$sry.treatment.p.adj <- NA
liver_T$chr.sry.treatment.p.adj <- NA

liver_T[,41:47] <- apply(liver_T[,41:47], 2, function(x){as.numeric(as.character(x))})
liver_E[,41:47] <- apply(liver_E[,41:47], 2, function(x){as.numeric(as.character(x))})

for(i in 41:47){
  liver_E[1:92,i+7] <- p.adjust(liver_E[1:92,i], method = 'fdr')
}
for(i in 41:47){
  liver_T[1:92,i+7] <- p.adjust(liver_T[1:92,i], method = 'fdr')
}
liver_T[90:96, 40:54]
write.csv(liver_T, "liver/network_liver_3way_T.csv", row.names = T, col.names = T)
write.csv(liver_E, "liver/network_liver_3way_E.csv", row.names = T, col.names = T)

##==========2-way===========
ing <- as.data.frame(t(ing))

ing_B <- ing[,which(ing["treatment",] == "B")]
ing_T <- ing[,which(ing["treatment",] == "T")]
ing_E <- ing[,which(ing["treatment",] == "E")]

ing_B$chr.p <- NA
ing_B$sry.p <- NA
ing_B$chr.sry.p <- NA

ing_T$chr.p <- NA
ing_T$sry.p <- NA
ing_T$chr.sry.p <- NA

ing_E$chr.p <- NA
ing_E$sry.p <- NA
ing_E$chr.sry.p <- NA

ing_B <- as.data.frame(t(ing_B))
ing_T <- as.data.frame(t(ing_T))
ing_E <- as.data.frame(t(ing_E))

ing_B$chr <- factor(ing_B$chr)
ing_B$sry <- factor(ing_B$sry)
ing_T$sry <- factor(ing_T$sry)
ing_T$chr <- factor(ing_T$chr)
ing_E$chr <- factor(ing_E$chr)
ing_E$sry <- factor(ing_E$sry)



for(i in 1:128){
  fm <- aov(ing_B[1:17,i] ~ chr * sry, data = ing_B[1:17,])
  ing_B[18:20, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:128){
  fm <- aov(unlist(ing_T[1:19,i]) ~ chr * sry, data = ing_T[1:19,])
  ing_T[20:22, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:128){
  fm <- aov(unlist(ing_E[1:19,i]) ~ chr * sry, data = ing_E[1:19,])
  ing_E[20:22, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

ing_B <- as.data.frame(t(ing_B))
ing_T <- as.data.frame(t(ing_T))
ing_E <- as.data.frame(t(ing_E))

ing_B$chr.p.adj <- NA
ing_B$sry.p.adj <- NA
ing_B$chr.sry.p.adj <- NA

ing_T$chr.p.adj <- NA
ing_T$sry.p.adj <- NA
ing_T$chr.sry.p.adj <- NA

ing_E$chr.p.adj <- NA
ing_E$sry.p.adj <- NA
ing_E$chr.sry.p.adj <- NA

ing_B$chr.p.adj[1:128] <- p.adjust(as.numeric(ing_B$chr.p[1:128]), method = 'fdr')
ing_B$sry.p.adj[1:128] <- p.adjust(as.numeric(ing_B$sry.p[1:128]), method = 'fdr')
ing_B$chr.sry.p.adj[1:128] <- p.adjust(as.numeric(ing_B$chr.sry.p[1:128]), method = 'fdr')

ing_T$chr.p.adj[1:128] <- p.adjust(as.numeric(ing_T$chr.p[1:128]), method = 'fdr')
ing_T$sry.p.adj[1:128] <- p.adjust(as.numeric(ing_T$sry.p[1:128]), method = 'fdr')
ing_T$chr.sry.p.adj[1:128] <- p.adjust(as.numeric(ing_T$chr.sry.p[1:128]), method = 'fdr')

ing_E$chr.p.adj[1:128] <- p.adjust(as.numeric(ing_E$chr.p[1:128]), method = 'fdr')
ing_E$sry.p.adj[1:128] <- p.adjust(as.numeric(ing_E$sry.p[1:128]), method = 'fdr')
ing_E$chr.sry.p.adj[1:128] <- p.adjust(as.numeric(ing_E$chr.sry.p[1:128]), method = 'fdr')

write.csv(ing_B, "ING/network_2wayANOVA_B.csv")
write.csv(ing_E, "ING/network_2wayANOVA_E.csv")
write.csv(ing_T, "ING/network_2wayANOVA_T.csv")

liver_B <- liver[which(liver$treatment == "B"),]
liver_T <- liver[which(liver$treatment == "T"),]
liver_E <- liver[which(liver$treatment == "E"),]

liver_B <- as.data.frame(t(liver_B))
liver_T <- as.data.frame(t(liver_T))
liver_E <- as.data.frame(t(liver_E))

liver_B$chr.p <- NA
liver_B$sry.p <- NA
liver_B$chr.sry.p <- NA

liver_T$chr.p <- NA
liver_T$sry.p <- NA
liver_T$chr.sry.p <- NA

liver_E$chr.p <- NA
liver_E$sry.p <- NA
liver_E$chr.sry.p <- NA

liver_B <- as.data.frame(t(liver_B))
liver_T <- as.data.frame(t(liver_T))
liver_E <- as.data.frame(t(liver_E))

liver_B$chr <- factor(liver_B$chr)
liver_B$sry <- factor(liver_B$sry)
liver_T$sry <- factor(liver_T$sry)
liver_T$chr <- factor(liver_T$chr)
liver_E$chr <- factor(liver_E$chr)
liver_E$sry <- factor(liver_E$sry)



for(i in 1:92){
  fm <- aov(liver_B[1:20,i] ~ chr * sry, data = liver_B[1:20,])
  liver_B[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:92){
  fm <- aov(liver_T[1:20,i] ~ chr * sry, data = liver_T[1:20,])
  liver_T[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}
#liver_E <- read.csv('liver/expr.lumi.gene(60_7466)_log2_0.01_dr0.3_2wayANOVA_E.csv')

for(i in 1:92){
  fm <- aov(liver_E[1:20,i] ~ chr * sry, data = liver_E[1:20,])
  liver_E[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

liver_B <- as.data.frame(t(liver_B))
liver_T <- as.data.frame(t(liver_T))
liver_E <- as.data.frame(t(liver_E))

liver_B$chr.p.adj <- NA
liver_B$sry.p.adj <- NA
liver_B$chr.sry.p.adj <- NA

liver_E$chr.p.adj <- NA
liver_E$sry.p.adj <- NA
liver_E$chr.sry.p.adj <- NA

liver_T$chr.p.adj <- NA
liver_T$sry.p.adj <- NA
liver_T$chr.sry.p.adj <- NA

# liver_B[,21:26] <- apply(liver_B[,21:26], 2, function(x){as.numeric(as.character(x))})
# liver_T[,21:26] <- apply(liver_T[,21:26], 2, function(x){as.numeric(as.character(x))})
# liver_E[,21:26] <- apply(liver_E[,21:26], 2, function(x){as.numeric(as.character(x))})

liver_B$chr.p.adj[1:92] <- p.adjust(liver_B$chr.p[1:92], method = 'fdr')
liver_B$sry.p.adj[1:92] <- p.adjust(liver_B$sry.p[1:92], method = 'fdr')
liver_B$chr.sry.p.adj[1:92] <- p.adjust(liver_B$chr.sry.p[1:92], method = 'fdr')

liver_T$chr.p.adj[1:92] <- p.adjust(liver_T$chr.p[1:92], method = 'fdr')
liver_T$sry.p.adj[1:92] <- p.adjust(liver_T$sry.p[1:92], method = 'fdr')
liver_T$chr.sry.p.adj[1:92] <- p.adjust(liver_T$chr.sry.p[1:92], method = 'fdr')

liver_E$chr.p.adj[1:92] <- p.adjust(liver_E$chr.p[1:92], method = 'fdr')
liver_E$sry.p.adj[1:92] <- p.adjust(liver_E$sry.p[1:92], method = 'fdr')
liver_E$chr.sry.p.adj[1:92] <- p.adjust(liver_E$chr.sry.p[1:92], method = 'fdr')

write.csv(liver_B, "liver/network_2wayANOVA_B.csv")
write.csv(liver_T, "liver/network_2wayANOVA_T.csv")
write.csv(liver_E, "liver/network_2wayANOVA_E.csv")

#==========pair-wise==========
ing <- t(ing)
ing <- as.data.frame(ing)

group <- as.data.frame(ing[1:55,129:132])
library(plyr)
group[35:40,]
group$chr <- revalue(group$chr, c("XY-" = "XY"))
group$sry <- revalue(group$sry, c("-" = "NULL"))
chr_g <- unique(group$chr)
sry_g <- unique(group$sry)
treatment_g <- unique(group$treatment)

ing <- t(ing[1:55, 1:128])
ing <- as.data.frame(ing)

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
                         file = paste0("ING/network_T-E_0.3_",chr,"_",sry,".xlsx"))
    
    tt2 <- topTable(raw.fit, coef = 2, n = Inf, adjust.method = "BH")
    tt2 <- tt2[order(tt2$adj.P.Val, decreasing = F),]
    tt2$gene <- rownames(tt2)
    openxlsx::write.xlsx(tt2, 
                         file = paste0("ING/network_E-B_0.3_",chr,"_",sry,".xlsx"))
    
    tt3 <- topTable(raw.fit, coef = 3, n = Inf, adjust.method = "BH")
    tt3 <- tt3[order(tt3$adj.P.Val, decreasing = F),]
    tt3$gene <- rownames(tt3)
    openxlsx::write.xlsx(tt3, 
                         file = paste0("ING/network_T-B_0.3_",chr,"_",sry,".xlsx")) 
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
                         file = paste0("ING/network_Sry-Null_0.3_",chr,"_",treatment,".xlsx"))
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
                         file = paste0("ING/network_XX-XY_0.3_new",sry,"_",treatment,".xlsx"))
    }
}

group <- as.data.frame(liver[1:60,93:96])
library(plyr)
group[35:40,]
group$chr <- revalue(group$chr, c("XY-" = "XY"))
group$sry <- revalue(group$sry, c("-" = "NULL"))
chr_g <- unique(group$chr)
sry_g <- unique(group$sry)
treatment_g <- unique(group$treatment)

liver <- t(liver[1:60, 1:92])
liver <- as.data.frame(liver)

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
                         file = paste0("liver/network_T-E_0.3_",chr,"_",sry,".xlsx"))
    
    tt2 <- topTable(raw.fit, coef = 2, n = Inf, adjust.method = "BH")
    tt2 <- tt2[order(tt2$adj.P.Val, decreasing = F),]
    tt2$gene <- rownames(tt2)
    openxlsx::write.xlsx(tt2, 
                         file = paste0("liver/network_E-B_0.3_",chr,"_",sry,".xlsx"))

    tt3 <- topTable(raw.fit, coef = 3, n = Inf, adjust.method = "BH")
    tt3 <- tt3[order(tt3$adj.P.Val, decreasing = F),]
    tt3$gene <- rownames(tt3)
    openxlsx::write.xlsx(tt3, 
                         file = paste0("liver/network_T-B_0.3_",chr,"_",sry,".xlsx"))
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
                         file = paste0("liver/network_Sry-Null_0.3_",chr,"_",treatment,".xlsx"))
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
                         file = paste0("liver/network_XX-XY_0.3_",sry,"_",treatment,".xlsx"))
  }
}