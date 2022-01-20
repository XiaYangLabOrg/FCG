ing <- read.table("ING_expr.lumi.gene(55_9400)_vst_0.01_dr0.3.txt", header = T)
ing <- t(ing)
ing <- as.data.frame(ing)
sample.info <- as.data.frame(read.table("Group_info_ING.txt", header = T))
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

ing_T <- ing[which(ing$treatment != "E"),]
ing_E <- ing[which(ing$treatment != "T"),]

library(tidyr)
ing_T <- split(ing_T, list(ing_T$chr, ing_T$sry))
ing_E <- split(ing_E, list(ing_E$chr, ing_E$sry))

for (j in 1:4) {
  for(i in 1:9400){
    nrow <- nrow(ing_E[[j]])
    fm <- aov(ing_E[[j]][,i] ~ treatment, data = ing_E[[j]])
    if(i == 1){
      ing_E[[j]][nrow + 1, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }else{
      ing_E[[j]][nrow, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }
    
  }
}

for (j in 1:4) {
  nrow <- nrow(ing_E[[j]])
  ing_E[[j]][nrow + 1, 1:9400] <- p.adjust(as.numeric(ing_E[[j]][nrow, 1:9400]), method = 'fdr')
  ing_E[[j]] <- as.data.frame(t(ing_E[[j]]))
}

for(j in 1:4){
  write.csv(ing_E[[j]], paste0("ing_E_", names(ing_E)[j],"1wayANOVA.csv"))
}

for (j in 1:4) {
  for(i in 1:9400){
    nrow <- nrow(ing_T[[j]])
    fm <- aov(ing_T[[j]][,i] ~ treatment, data = ing_T[[j]])
    if(i == 1){
      ing_T[[j]][nrow + 1, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }else{
      ing_T[[j]][nrow, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }
    
  }
}

for (j in 1:4) {
  nrow <- nrow(ing_T[[j]])
  ing_T[[j]][nrow + 1, 1:9400] <- p.adjust(as.numeric(ing_T[[j]][nrow, 1:9400]), method = 'fdr')
  ing_T[[j]] <- as.data.frame(t(ing_T[[j]]))
}

for(j in 1:4){
  write.csv(ing_T[[j]], paste0("ing_T_", names(ing_T)[j],"1wayANOVA.csv"))
}

#==============liver==========
liver <- read.table("liver_expr.lumi.gene(60_7466)_vst_0.01_dr0.3.txt", header = T)
# liver <- as.data.frame(liver)
liver <- t(liver)
liver <- as.data.frame(liver)
sample.info <- as.data.frame(read.table("Group_info_liver.txt", header = T))
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

liver_T <- liver[which(liver$treatment != "E"),]
liver_E <- liver[which(liver$treatment != "T"),]

library(tidyr)
liver_T <- split(liver_T, list(liver_T$chr, liver_T$sry))
liver_E <- split(liver_E, list(liver_E$chr, liver_E$sry))

for (j in 1:4) {
  for(i in 1:7466){
    nrow <- nrow(liver_E[[j]])
    fm <- aov(liver_E[[j]][,i] ~ treatment, data = liver_E[[j]])
    if(i == 1){
      liver_E[[j]][nrow + 1, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }else{
      liver_E[[j]][nrow, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }
    
  }
}

for (j in 1:4) {
  nrow <- nrow(liver_E[[j]])
  liver_E[[j]][nrow + 1, 1:7466] <- p.adjust(as.numeric(liver_E[[j]][nrow, 1:7466]), method = 'fdr')
  liver_E[[j]] <- as.data.frame(t(liver_E[[j]]))
}

for(j in 1:4){
  write.csv(liver_E[[j]], paste0("liver_E_", names(liver_E)[j],"1wayANOVA.csv"))
}

for (j in 1:4) {
  for(i in 1:7466){
    nrow <- nrow(liver_T[[j]])
    fm <- aov(liver_T[[j]][,i] ~ treatment, data = liver_T[[j]])
    if(i == 1){
      liver_T[[j]][nrow + 1, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }else{
      liver_T[[j]][nrow, i] <- summary(fm)[[1]]$`Pr(>F)`[1]
    }
    
  }
}

for (j in 1:4) {
  nrow <- nrow(liver_T[[j]])
  liver_T[[j]][nrow + 1, 1:7466] <- p.adjust(as.numeric(liver_T[[j]][nrow, 1:7466]), method = 'fdr')
  liver_T[[j]] <- as.data.frame(t(liver_T[[j]]))
}

for(j in 1:4){
  write.csv(liver_T[[j]], paste0("liver_T_", names(liver_T)[j],"1wayANOVA.csv"))
}