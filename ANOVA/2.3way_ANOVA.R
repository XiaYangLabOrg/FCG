setwd("G:/UCLA/FCG/")
#=========ING 3way ANOVA=========
ing <- read.table("ING/ING_expr.lumi.gene(55_9400)_log2_0.01_dr0.3.txt", header = T)
ing <- t(ing)
ing <- as.data.frame(ing)
sample.info <- as.data.frame(read.table("ING/Group_info_ING.txt", header = T))
# sample.info$Sample.ID <- paste0("X", sample.info$Sample.ID)

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

for(i in 1:9400){
  fm <- aov(unlist(ing_T[1:36,i]) ~ chr * sry * treatment, data = ing_T[1:36,])
  ing_T[37:43, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}

for(i in 1:9400){
  fm <- aov(unlist(ing_E[1:36,i]) ~ chr * sry * treatment, data = ing_E[1:36,])
  ing_E[37:43, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}


#==========liver 3way ANOVA=========
liver <- read.table("liver/liver_expr.lumi.gene(60_7466)_log2_0.01_dr0.3.txt", header = T)
# liver <- as.data.frame(liver)
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

for(i in 1:7466){
  fm <- aov(unlist(liver_T[1:40,i]) ~ chr * sry * treatment, data = liver_T[1:40,])
  liver_T[41:47, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}
for(i in 1:7466){
  fm <- aov(unlist(liver_E[1:40,i]) ~ chr * sry * treatment, data = liver_E[1:40,])
  liver_E[41:47, i] <- summary(fm)[[1]]$`Pr(>F)`[1:7]
}


#======ING adjust p value====
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
  ing_E[1:9400,i+7] <- p.adjust(ing_E[1:9400,i], method = 'fdr')
}
for(i in 37:43){
  ing_T[1:9400,i+7] <- p.adjust(ing_T[1:9400,i], method = 'fdr')
}
ing_E[8900:8907, 37:50]

write.csv(ing_T, "ING/expr.lumi.gene(55_9400)_vst_0.01_dr0.45_3wayANOVA_T.csv", row.names = T, col.names = T)
write.csv(ing_E, "ING/expr.lumi.gene(55_9400)_vst_0.01_dr0.45_3wayANOVA_E.csv", row.names = T, col.names = T)

#======liver adjust p value=====
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
  liver_E[1:7466,i+7] <- p.adjust(liver_E[1:7466,i], method = 'fdr')
}
for(i in 41:47){
  liver_T[1:7466,i+7] <- p.adjust(liver_T[1:7466,i], method = 'fdr')
}
liver_T[7095:7101, 40:54]

write.csv(liver_T, "liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.3_3wayANOVA_T.csv", row.names = T, col.names = T)
write.csv(liver_E, "liver/expr.lumi.gene(60_7466)_vst_0.01_dr0.3_3wayANOVA_E.csv", row.names = T, col.names = T)
