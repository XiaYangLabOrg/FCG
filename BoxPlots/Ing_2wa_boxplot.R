library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggeasy)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)


#####====== T- Sry:Chr ===========

# No significant DEGs 


#####====== E- Sry:Chr ========
Ing_2wa_E_gonad_chr_DEGs <- read_excel("Data/Ing_2wa_E_gonad_chr_DEGs.xlsx")
g.c.E.reversed <- t(Ing_2wa_E_gonad_chr_DEGs)
g.c.E.reversed <- g.c.E.reversed[-c(21:26),]
rownames(g.c.E.reversed)[1] <- "Samples"
colnames(g.c.E.reversed) <- g.c.E.reversed[1,]
g.c.E.reversed <- g.c.E.reversed[-1,]
g.c.E.samples <- rownames(g.c.E.reversed)
g.c.E.samples <- as.data.frame(g.c.E.samples)
rownames(g.c.E.samples) <- g.c.E.samples[,1]
g.c.E.reversed <- cbind(g.c.E.reversed,g.c.E.samples)
g.c.E.reversed$genotype <- "temp"
rownames(g.c.E.reversed) <- c(1:19)
g.c.E.reversed$genotype[1:5] <- "XYM"
g.c.E.reversed$genotype[6:9] <- "XXM"
g.c.E.reversed$genotype[10:14] <- "XYF"
g.c.E.reversed$genotype[14:19] <- "XXF"

########### Map7d3
# Map7d3.Egc <- g.c.E.reversed
# colnames(Map7d3.Egc)[1] <- "Expression_level"
# Map7d3.Egc$genotype <- as.factor(Map7d3.Egc$genotype)
# Map7d3.Egc$Expression_level <- as.numeric(Map7d3.Egc$Expression_level)
# png(file = "Plots/Ing_2wa_E_gc_Map7d3.png")
# Map7d3.Egc_plot <- ggplot(Map7d3.Egc,aes(x = genotype, y = Expression_level))+ 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Map7d3 in Adipose 2WA (gonad:chr) [E] ", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Map7d3.Egc_plot
# dev.off()



Map7d3.Egc <- g.c.E.reversed[,c("Map7d3","genotype")]
colnames(Map7d3.Egc)[1] <- "Expression_level"
Map7d3.Egc$genotype <- as.factor(Map7d3.Egc$genotype)
Map7d3.Egc$Expression_level <- as.numeric(Map7d3.Egc$Expression_level)

Map7d3.Egc_plot <- ggplot(Map7d3.Egc, aes(x= genotype, y= Expression_level)) + 
  geom_boxplot(outlier.shape = NA, fill = "#d04290") + 
  labs(title = "Map7d3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4, fill = "#d04290") 

## pairwise comparisons
Map7d3.Egc.stat <- Map7d3.Egc %>% 
  t_test(Expression_level ~ genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Map7d3.Egc.stat

Map7d3.Egc.stat <- Map7d3.Egc.stat %>%
  add_xy_position(x = "genotype")

Map7d3.Egc_plot <- Map7d3.Egc_plot+
  stat_pvalue_manual(
    Map7d3.Egc.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_2wa_E_gc/Map7d3.png", Map7d3.Egc_plot, width = 5, height = 5, dpi = 300)

#####====== B- Sry:Chr ========
Ing_2wa_B_gonad_chr_DEGs <- read_excel("Data/Ing_2wa_B_gonad_chr_DEGs.xlsx")
g.c.B.reversed <- t(Ing_2wa_B_gonad_chr_DEGs)
rownames(g.c.B.reversed)[1] <- "Samples" 
colnames(g.c.B.reversed) <- g.c.B.reversed[1,]
g.c.B.reversed <- g.c.B.reversed[-1,]
g.c.B.samples <- rownames(g.c.B.reversed)
g.c.B.samples <- as.data.frame(g.c.B.samples)
rownames(g.c.B.samples) <- g.c.B.samples[,1]
g.c.B.reversed <- cbind(g.c.B.reversed,g.c.B.samples)
g.c.B.reversed$genotype <- "temp"
rownames(g.c.B.reversed) <- c(1:17)
g.c.B.reversed <- g.c.B.reversed[order(g.c.B.reversed$group),]
g.c.B.reversed$genotype[1:3] <- "XXF"
g.c.B.reversed$genotype[4:8] <- "XXM"
g.c.B.reversed$genotype[9:13] <- "XYF"
g.c.B.reversed$genotype[13:17] <- "XYM"

# Mpp5

Mpp5.Bgc <- g.c.B.reversed[,c("Mpp5","genotype")]
colnames(Mpp5.Bgc)[1] <- "Expression_level"
Mpp5.Bgc$genotype <- as.factor(Mpp5.Bgc$genotype)
Mpp5.Bgc$Expression_level <- as.numeric(Mpp5.Bgc$Expression_level)

Mpp5.Bgc_plot <- ggplot(Mpp5.Bgc, aes(x= genotype, y= Expression_level)) + 
  geom_boxplot(outlier.shape = NA, fill = "white") + 
  labs(title = "Mpp5", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4, fill = "white") 

## pairwise comparisons
Mpp5.Bgc.stat <- Mpp5.Bgc %>% 
  t_test(Expression_level ~ genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Mpp5.Bgc.stat

Mpp5.Bgc.stat <- Mpp5.Bgc.stat %>%
  add_xy_position(x = "genotype")

Mpp5.Bgc_plot <- Mpp5.Bgc_plot+
  stat_pvalue_manual(
    Mpp5.Bgc.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_2wa_B_gc/Mpp5.png", Mpp5.Bgc_plot, width = 5, height = 5, dpi = 300)

# Sh3d4
Sh3d4.Bgc <- g.c.B.reversed[,c("Sh3d4","genotype")]
colnames(Sh3d4.Bgc)[1] <- "Expression_level"
Sh3d4.Bgc$genotype <- as.factor(Sh3d4.Bgc$genotype)
Sh3d4.Bgc$Expression_level <- as.numeric(Sh3d4.Bgc$Expression_level)

Sh3d4.Bgc_plot <- ggplot(Sh3d4.Bgc, aes(x= genotype, y= Expression_level)) + 
  geom_boxplot(outlier.shape = NA, fill = "white") + 
  labs(title = "Sh3d4", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4, fill = "white") 

## pairwise comparisons
Sh3d4.Bgc.stat <- Sh3d4.Bgc %>% 
  t_test(Expression_level ~ genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Sh3d4.Bgc.stat

Sh3d4.Bgc.stat <- Sh3d4.Bgc.stat %>%
  add_xy_position(x = "genotype")

Sh3d4.Bgc_plot <- Sh3d4.Bgc_plot+
  stat_pvalue_manual(
    Sh3d4.Bgc.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_2wa_B_gc/Sh3d4.png", Sh3d4.Bgc_plot, width = 5, height = 5, dpi = 300)


##### arrange

ing.2wa.B.gc.plots <- ggarrange(Mpp5.Bgc_plot, Sh3d4.Bgc_plot,
                                labels = c("A", "B"),
                                ncol = 2, nrow = 1)
ggsave("Plots_new/Ing_2wa_B_gc/Combined.png", ing.2wa.B.gc.plots,width = 10, height = 5, dpi = 600)
