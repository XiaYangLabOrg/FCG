library(limma)
library(lumi)
setwd("G:/UCLA/FCG/")
targets.ing <- read.table("ING/Group_info_ING.txt", header = T)

##=========using lumi===========

##===========ING==========
ing.lumi <- lumiR.batch("ING/Sample_Probe_Profile.txt", convertNuID = F)
ing.lumi <- addControlData2lumi("ING/Control_Probe_Profile.txt", ing.lumi)
ing.lumi@assayData$exprs[1:5,1:5]

sampleNames(ing.lumi) <- pData(ing.lumi)$sampleID
badsamp <- c("5613711027_E", "5613711014_A", "5613711014_G", "5719882029_E", "5719882029_H")
toremove <- pData(ing.lumi)$sampleID %in% badsamp
ing.lumi <- ing.lumi[, -which(toremove)]

summary(ing.lumi, "QC")
plot(ing.lumi, what = 'density', legend = F)
plotCDF(ing.lumi.NQ, reverse = TRUE, addLegend = F)
plot(ing.lumi.NQ, what = 'sampleRelation')

ing.lumi.NQ <- lumiExpresso(ing.lumi, bg.correct = T, bgcorrect.param = list(method = 'bgAdjust'),
                            variance.stabilize = T, varianceStabilize.param = list(), normalize = T,
                            normalize.param = list(method = 'quantile'), QC.evaluation = T, QC.param = list(),
                            verbose = T)
summary(ing.lumi.NQ, 'QC')
ing.lumi.NQ@assayData$exprs[1:5,1:5]

ing.lumi.expr <- exprs(ing.lumi.NQ)
dim(ing.lumi.expr)
ing.lumi.p <- ing.lumi.NQ@assayData$detection
featureData <- fData(ing.lumi.NQ)
detect.count <- detectionCall(ing.lumi.NQ, Th = 0.01, type = "probe")
detect.rate <- detect.count/ncol(ing.lumi.expr)
detectP <- detectionCall(ing.lumi.NQ, Th = 0.01, type = "matrix")
detectP[detectP == "P"] <- T
detectP[detectP == "A"] <- F
mode(detectP) <- "logical"
ing.lumi.expr[!detectP] <- 0
expr.lumi <- ing.lumi.expr[detect.rate >= 0.30,]

boxplot(expr.lumi,range = 0)
rownames(expr.lumi) <- ing.lumi@featureData@data$SYMBOL[match(rownames(expr.lumi), ing.lumi@featureData@data$PROBE_ID)]

expr.lumi <- apply(expr.lumi, 2, function(x){tapply(x, row.names(expr.lumi), mean)})

expr.lumi <- as.data.frame(expr.lumi)

expr.lumi.gg <- as.data.frame(t(expr.lumi))
expr.lumi.gg$ID <- rownames(expr.lumi.gg)
expr.lumi.gg$group <- targets.ing$Sample.Group[match(expr.lumi.gg$ID, targets.ing$Sample.ID)]
expr.lumi.gg$chr[expr.lumi.gg$group %in% c("XX_B", "XX_E", "XX_T")] <- "XX"
expr.lumi.gg$chr[expr.lumi.gg$group %in% c("XXSry_B", "XXSry_E", "XXSry_T")] <- "XXSry"
expr.lumi.gg$chr[expr.lumi.gg$group %in% c("XY-_B", "XY-_E", "XY-_T")] <- "XY-"
expr.lumi.gg$chr[expr.lumi.gg$group %in% c("XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY-Sry"

write.table(expr.lumi, "ING_expr(55_9400)_vst_0.01_dr0.5.txt", row.names = T, col.names = T)

boxplot(expr.lumi, main = "boxplot of lumi generated data", range = 0)
expr.lumi.gg$ID <- factor(expr.lumi.gg$ID, levels = expr.lumi.gg$ID, ordered = T)

ggplot(expr.lumi.gg[order(expr.lumi.gg$group),], aes(x = ID, y = Eif2s3x, fill = group))+
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8))

ggplot(expr.lumi, aes(x = expr.lumi$`5719882029_E`, y = expr.lumi$`5613711027_B`))+
  geom_point()

barplot(expr.lumi["Jarid1d",], col = as.factor(targets.ing$Sample.Group), 
        main = "Jarid1d - using lumi for norm&bkgd", names.arg = targets.ing$Sample.Group)

#tsne
set.seed(1)
library(Rtsne)
tsne <- Rtsne(t(expr.lumi), dims = 2, perplexity = 7)
tsne_for_plot <- as.data.frame(tsne$Y)
tsne_for_plot$ID <- colnames(expr.lumi)
tsne_for_plot$group.original <- targets.ing$Sample.Group[match(colnames(expr.lumi), targets.ing$Sample.ID)]
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
tsne_for_plot$chromosome[tsne_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
tsne_for_plot$chromosome[tsne_for_plot$group.original %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY"
tsne_for_plot$sex[tsne_for_plot$group.original %in% c("XXSry_B", "XXSry_E", "XXSry_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "M"
tsne_for_plot$sex[tsne_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "F"
library(ggplot2)
assumed_outlier <- tsne_for_plot$ID %in% c("5613711027_E", "5613711014_A", "5613711014_G")
ggplot(tsne_for_plot, aes(x = V1, y = V2, col = chromosome)) +
  geom_point(size = 5)# + #+
  #, shape = factor(tsne_for_plot$group.original)) #+
  #geom_text(label = tsne_for_plot$ID)
  # xlim(-7,7) +
  # ylim(-10,9)

#pca
pca <- prcomp(t(expr.lumi), center = T)
plot(pca, type = "l")
summary(pca)
pca_for_plot <- as.data.frame(pca$x)
pca_for_plot$ID <- colnames(expr.lumi)
pca_for_plot$group.original <- targets.ing$Sample.Group[match(colnames(expr.lumi), targets.ing$Sample.ID)]
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
pca_for_plot$chromosome[pca_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
pca_for_plot$chromosome[pca_for_plot$group.original %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY"
pca_for_plot$sex[pca_for_plot$group.original %in% c("XXSry_B", "XXSry_E", "XXSry_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "M"
pca_for_plot$sex[pca_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "F"
ggplot(pca_for_plot, aes(x = PC1, y = PC2, col = treatment))+
  geom_point(size = 4) #+
  #geom_text(label = pca_for_plot$ID)

#=============liver=========
setwd("G:/UCLA/FCG/liver/")
targets.liver <- read.table("Group_info_liver.txt", header = T)

liver.lumi <- lumiR.batch("Sample_Probe_Profile.txt", convertNuID = F)
liver.lumi <- addControlData2lumi("Control_Probe_Profile.txt", liver.lumi)
liver.lumi@assayData$exprs[1:5,1:5]
summary(liver.lumi, "QC")
plot(liver.lumi.NQ, what = 'density', legend = F)
plotCDF(liver.lumi.NQ, reverse = TRUE, addLegend = F)
plot(liver.lumi, what = 'sampleRelation')

liver.lumi.NQ <- lumiExpresso(liver.lumi, bg.correct = T, bgcorrect.param = list(method = 'bgAdjust'),
                            variance.stabilize = T, varianceStabilize.param = list(), normalize = T,
                            normalize.param = list(method = 'quantile'), QC.evaluation = T, QC.param = list(),
                            verbose = T)#
summary(liver.lumi.NQ, 'QC')
liver.lumi.NQ@assayData$exprs[1:5,1:5]

liver.lumi.expr <- exprs(liver.lumi.NQ)
dim(liver.lumi.expr)
liver.lumi.p <- liver.lumi.NQ@assayData$detection
featureData <- fData(liver.lumi.NQ)
detect.count <- detectionCall(liver.lumi.NQ, Th = 0.01, type = "probe")
detect.rate <- detect.count/ncol(liver.lumi.expr)
detectP <- detectionCall(liver.lumi.NQ, Th = 0.01, type = "matrix")
detectP[detectP == "P"] <- T
detectP[detectP == "A"] <- F
mode(detectP) <- "logical"
liver.lumi.expr[!detectP] <- 0
expr.lumi.liver <- liver.lumi.expr[detect.rate >= 0.30,]
boxplot(expr.lumi.liver, range = 0)
rownames(expr.lumi.liver) <- liver.lumi@featureData@data$SYMBOL[match(rownames(expr.lumi.liver), liver.lumi@featureData@data$PROBE_ID)]

expr.lumi.liver <- apply(expr.lumi.liver, 2, function(x){tapply(x, row.names(expr.lumi.liver), mean)})
expr.lumi.liver <- as.data.frame(expr.lumi.liver)

write.table(expr.lumi.liver, "liver_expr(60_7466)_vst_log2_0.01_dr0.3.txt", row.names = T, col.names = T)

expr.lumi.liver.gg <- as.data.frame(t(expr.lumi.liver))
expr.lumi.liver.gg$ID <- rownames(expr.lumi.liver.gg)
expr.lumi.liver.gg$group <- targets.liver$Sample.Group[match(expr.lumi.liver.gg$ID, targets.liver$Sample.ID)]
expr.lumi.liver.gg$chr[expr.lumi.liver.gg$group %in% c("XX_B", "XX_E", "XX_T")] <- "XX"
expr.lumi.liver.gg$chr[expr.lumi.liver.gg$group %in% c("XXSry_B", "XXSry_E", "XXSry_T")] <- "XXSry"
expr.lumi.liver.gg$chr[expr.lumi.liver.gg$group %in% c("XY-_B", "XY-_E", "XY-_T")] <- "XY-"
expr.lumi.liver.gg$chr[expr.lumi.liver.gg$group %in% c("XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY-Sry"

expr.lumi.liver.gg$ID <- factor(expr.lumi.liver.gg$ID, levels = expr.lumi.liver.gg$ID, ordered = T)

ggplot(expr.lumi.liver.gg[expr.lumi.liver.gg$group,], aes(x = ID, y = Eif2s3x, fill = group))+
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8))

barplot(expr.lumi.liver["Jarid1d",], col = as.factor(targets.liver$Sample.Group), main = "Jarid1d - using lumi for norm&bkgd", names.arg = targets.liver$Sample.Group)

set.seed(2)
library(Rtsne)
tsne <- Rtsne(t(expr.lumi.liver), dims = 2, perplexity = 7)
tsne_for_plot <- as.data.frame(tsne$Y)
tsne_for_plot$ID <- colnames(expr.lumi.liver)
tsne_for_plot$group.original <- targets.liver$Sample.Group[match(colnames(expr.lumi.liver), targets.liver$Sample.ID)]
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
tsne_for_plot$treatment[tsne_for_plot$group.original %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
tsne_for_plot$chromosome[tsne_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
tsne_for_plot$chromosome[tsne_for_plot$group.original %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY"
tsne_for_plot$sex[tsne_for_plot$group.original %in% c("XXSry_B", "XXSry_E", "XXSry_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "M"
tsne_for_plot$sex[tsne_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "F"
library(ggplot2)
assumed_outlier <- tsne_for_plot$ID %in% c("5719882028_B", "5613703028_E")
ggplot(tsne_for_plot, aes(x = V1, y = V2, col = sex)) +
  geom_point(size = 5)# +
  #ggtitle("perplexity = 7") #+
#, shape = factor(tsne_for_plot$group.original)) #+
#geom_text(label = tsne_for_plot$ID)
# xlim(-7,7) +
# ylim(-10,9)


#pca
pca <- prcomp(t(expr.lumi.liver), center = T)
plot(pca, type = "l")
summary(pca)
pca_for_plot <- as.data.frame(pca$x)
pca_for_plot$ID <- colnames(expr.lumi.liver)
pca_for_plot$group.original <- targets.liver$Sample.Group[match(colnames(expr.lumi.liver), targets.liver$Sample.ID)]
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_B", "XXSry_B", "XY-_B", "XY-Sry_B")] <- "B"
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_E", "XXSry_E", "XY-_E", "XY-Sry_E")] <- "E"
pca_for_plot$treatment[pca_for_plot$group.original %in% c("XX_T", "XXSry_T", "XY-_T", "XY-Sry_T")] <- "T"
pca_for_plot$chromosome[pca_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XXSry_B", "XXSry_E", "XXSry_T")] <- "XX"
pca_for_plot$chromosome[pca_for_plot$group.original %in% c("XY-_B", "XY-_E", "XY-_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "XY"
pca_for_plot$sex[pca_for_plot$group.original %in% c("XXSry_B", "XXSry_E", "XXSry_T", "XY-Sry_B", "XY-Sry_E", "XY-Sry_T")] <- "M"
pca_for_plot$sex[pca_for_plot$group.original %in% c("XX_B", "XX_E", "XX_T", "XY-_B", "XY-_E", "XY-_T")] <- "F"
ggplot(pca_for_plot, aes(x = PC1, y = PC2, col = treatment))+
  geom_point(size = 4) #+
  #geom_text(label = pca_for_plot$ID)


