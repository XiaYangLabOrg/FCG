heatmap <- read.csv("module_for heatmap.csv", stringsAsFactors = F,header = F)
library(ComplexHeatmap)
library(circlize)
num <- heatmap[3:94,4:46]
num <- apply(num, 2, function(x){as.numeric(x)})
rownames(num) <- heatmap[3:94,1]
colnames(num) <- heatmap[2,4:46]

annotation <- data.frame(test = c(rep("3-way ANOVA",14), rep("2-way ANOVA",9),rep("pair-wise", 20)),
                         background = as.character(heatmap[1,4:46]))
annotation$background <- sub(x = annotation$background,pattern = "3wayANOVA_liver_",replacement = "")
annotation$background <- sub(x = annotation$background,pattern = "2wayANOVA_",replacement = "")

GO <- heatmap[3:94,2]
GO <- GO[-which(apply(num, 1, min)>0.05)]

num <- num[-which(apply(num, 1, min)>0.05),]
num_log <- -log(num)

colnames(num_log) <- sub('treatment','horm',colnames(num_log))

background <- HeatmapAnnotation(background = anno_text(annotation[,2], rot = 90, just = "left",offset = unit(1, "npc") - unit(2, "mm")),show_annotation_name = TRUE)
test <- HeatmapAnnotation(df = annotation[,1,drop = F],show_annotation_name = TRUE,
                          factor = anno_text(colnames(num_log), rot = 90, just = "right",offset = unit(1, "npc") - unit(2, "mm")))
# col = colorspace::rainbow_hcl)
# text = anno_text(annotation[,2], rot = 45, just = "left", offset = unit(2, "mm")))


hub <- heatmap[3:94,3]
hubd <- gsub("\\(..\\)","",hub)
heat <- Heatmap(num_log,cluster_rows = T,cluster_columns = T,col = colorRamp2(c(1.3, 5, 45), c("white","yellow",rgb(87.5, 5.5, 5.5,maxColorValue = 100))),name = "-log(P-Value)",
                top_annotation = background, bottom_annotation = test, color_space = "RGB", show_column_names = F,
                show_row_dend=F, show_column_dend = F)
ha <- rowAnnotation(text = row_anno_text(GO,just = 0))
draw(ha)
heat
pdf("heatmap_p_liver_sub.pdf",height = 11, width = 18)
draw(heat + ha, padding = unit(c(30,2,40,130),"mm"),heatmap_legend_side = "left",annotation_legend_side = "left")

dev.off()