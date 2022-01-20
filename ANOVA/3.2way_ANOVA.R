##===========ING==========
ing <- read.csv("ING/expr.lumi.gene(55_9400)_vst_0.01_dr0.3_3wayANOVA.csv", header = T, stringsAsFactors = F)

ing_B <- ing[,which(ing["treatment",] == "B")]
ing_T <- ing[,which(ing["treatment",] == "T")]
ing_E <- ing[,which(ing["treatment",] == "E")]

ing_B <- as.data.frame(t(ing_B))
ing_T <- as.data.frame(t(ing_T))
ing_E <- as.data.frame(t(ing_E))

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

#=======2way ANOVA=======
for(i in 1:9400){
  fm <- aov(ing_B[1:17,i] ~ chr * sry, data = ing_B[1:17,])
  ing_B[18:20, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:9400){
  fm <- aov(unlist(ing_T[1:19,i]) ~ chr * sry, data = ing_T[1:19,])
  ing_T[20:22, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:9400){
  fm <- aov(unlist(ing_E[1:19,i]) ~ chr * sry, data = ing_E[1:19,])
  ing_E[20:22, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

ing_B <- as.data.frame(t(ing_B))
ing_T <- as.data.frame(t(ing_T))
ing_E <- as.data.frame(t(ing_E))

#=======p value adjustment======
ing_B$chr.p.adj <- NA
ing_B$sry.p.adj <- NA
ing_B$chr.sry.p.adj <- NA

ing_T$chr.p.adj <- NA
ing_T$sry.p.adj <- NA
ing_T$chr.sry.p.adj <- NA

ing_E$chr.p.adj <- NA
ing_E$sry.p.adj <- NA
ing_E$chr.sry.p.adj <- NA

ing_B$chr.p.adj[1:9400] <- p.adjust(as.numeric(ing_B$chr.p[1:9400]), method = 'fdr')
ing_B$sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_B$sry.p[1:9400]), method = 'fdr')
ing_B$chr.sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_B$chr.sry.p[1:9400]), method = 'fdr')

ing_T$chr.p.adj[1:9400] <- p.adjust(as.numeric(ing_T$chr.p[1:9400]), method = 'fdr')
ing_T$sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_T$sry.p[1:9400]), method = 'fdr')
ing_T$chr.sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_T$chr.sry.p[1:9400]), method = 'fdr')

ing_E$chr.p.adj[1:9400] <- p.adjust(as.numeric(ing_E$chr.p[1:9400]), method = 'fdr')
ing_E$sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_E$sry.p[1:9400]), method = 'fdr')
ing_E$chr.sry.p.adj[1:9400] <- p.adjust(as.numeric(ing_E$chr.sry.p[1:9400]), method = 'fdr')

write.csv(ing_B, "ING/expr.lumi.gene(55_9400).vst.0.01.dr0.3_2wayANOVA_B.csv")
write.csv(ing_E, "ING/expr.lumi.gene(55_9400).vst.0.01.dr0.3_2wayANOVA_E.csv")
write.csv(ing_T, "ING/expr.lumi.gene(55_9400).vst.0.01.dr0.3_2wayANOVA_T.csv")

#=============liver===========
liver <- read.table("liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.3.txt")
sample.info <- as.data.frame(read.table("liver/Targets.new.txt", header = T))
sample.info$Sample.ID <- paste0("X", sample.info$Sample.ID)

liver <- as.data.frame(t(liver))

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

rownames(liver) <- liver$X
liver <- liver[,-1]
liver <- as.data.frame(t(liver))

liver_B <- liver[which(liver$treatment == "B"),]
liver_T <- liver[which(liver$treatment == "T"),]
liver_E <- liver[which(liver$treatment == "E"),]

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

#==========2way ANOVA=====
for(i in 1:7466){
  fm <- aov(liver_B[1:20,i] ~ chr * sry, data = liver_B[1:20,])
  liver_B[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

for(i in 1:7466){
  fm <- aov(liver_T[1:20,i] ~ chr * sry, data = liver_T[1:20,])
  liver_T[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}
#liver_E <- read.csv('liver/expr.lumi.gene(60_7466)_log2_0.01_dr0.3_2wayANOVA_E.csv')

for(i in 1:7466){
  fm <- aov(liver_E[1:20,i] ~ chr * sry, data = liver_E[1:20,])
  liver_E[21:23, i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}

#========p value adjustemtn=====
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

liver_B[,21:26] <- apply(liver_B[,21:26], 2, function(x){as.numeric(as.character(x))})
liver_T[,21:26] <- apply(liver_T[,21:26], 2, function(x){as.numeric(as.character(x))})
liver_E[,21:26] <- apply(liver_E[,21:26], 2, function(x){as.numeric(as.character(x))})

liver_B$chr.p.adj[1:7466] <- p.adjust(liver_B$chr.p[1:7466], method = 'fdr')
liver_B$sry.p.adj[1:7466] <- p.adjust(liver_B$sry.p[1:7466], method = 'fdr')
liver_B$chr.sry.p.adj[1:7466] <- p.adjust(liver_B$chr.sry.p[1:7466], method = 'fdr')

liver_T$chr.p.adj[1:7466] <- p.adjust(liver_T$chr.p[1:7466], method = 'fdr')
liver_T$sry.p.adj[1:7466] <- p.adjust(liver_T$sry.p[1:7466], method = 'fdr')
liver_T$chr.sry.p.adj[1:7466] <- p.adjust(liver_T$chr.sry.p[1:7466], method = 'fdr')

liver_E$chr.p.adj[1:7466] <- p.adjust(liver_E$chr.p[1:7466], method = 'fdr')
liver_E$sry.p.adj[1:7466] <- p.adjust(liver_E$sry.p[1:7466], method = 'fdr')
liver_E$chr.sry.p.adj[1:7466] <- p.adjust(liver_E$chr.sry.p[1:7466], method = 'fdr')

write.csv(liver_B, "liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.5_2wayANOVA_B.csv")
write.csv(liver_T, "liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.5_2wayANOVA_T.csv")
write.csv(liver_E, "liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.5_2wayANOVA_E.csv")