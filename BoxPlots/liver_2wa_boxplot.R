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


#####===== E- Sry:Chr ===========


Liver_2wa_E_gc_DEGs <- read_excel("Data/Liver_2wa_E_gonad_chr_DEGs.xlsx")
g.c.E.reversed <- t(Liver_2wa_E_gc_DEGs)
colnames(g.c.E.reversed) <- g.c.E.reversed[1,]
g.c.E.reversed <- g.c.E.reversed[-1,]
g.c.E.reversed <- as.data.frame(g.c.E.reversed)
g.c.E.reversed$genotype <- "temp"
rownames(g.c.E.reversed) <- c(1:26)
g.c.E.reversed <- g.c.E.reversed[-c(21:26),]
g.c.E.reversed$genotype[1:5] <- "XYM"
g.c.E.reversed$genotype[6:10] <- "XXM"
g.c.E.reversed$genotype[11:15] <- "XYF"
g.c.E.reversed$genotype[15:20] <- "XXF"

########### Plekhg6  ###############

# Plekhg6.Tgc <- g.c.E.reversed[c("Plekhg6","genotype")]
# colnames(Plekhg6.Tgc)[1] <- "Expression_level"
# Plekhg6.Tgc$genotype <- as.factor(Plekhg6.Tgc$genotype)
# Plekhg6.Tgc$Expression_level <- as.numeric(Plekhg6.Tgc$Expression_level)
# png(file = "Plots/Liver_2wa_E_gc_Plekhg6.png")
# Plekhg6.Tgc_plot <- ggplot(Plekhg6.Tgc,aes(x = genotype, y = Expression_level))+ 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Plekhg6 in Liver 2WA (gonad:chr) [E] ", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Plekhg6.Tgc_plot
# dev.off()


Plekhg6.Egc <- g.c.E.reversed[,c("Plekhg6","genotype")]
colnames(Plekhg6.Egc)[1] <- "Expression_level"
Plekhg6.Egc$genotype <- as.factor(Plekhg6.Egc$genotype)
Plekhg6.Egc$Expression_level <- as.numeric(Plekhg6.Egc$Expression_level)

Plekhg6.Egc_plot <- ggplot(Plekhg6.Egc, aes(x= genotype, y= Expression_level)) + 
  geom_boxplot(outlier.shape = NA, fill = "#d04290") + 
  labs(title = "Plekhg6", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4, fill = "#d04290") 

## pairwise comparisons
Plekhg6.Egc.stat <- Plekhg6.Egc %>% 
  t_test(Expression_level ~ genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Plekhg6.Egc.stat

Plekhg6.Egc.stat <- Plekhg6.Egc.stat %>%
  add_xy_position(x = "genotype")

Plekhg6.Egc_plot <- Plekhg6.Egc_plot+
  stat_pvalue_manual(
    Plekhg6.Egc.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Liver_2wa_E_gc/Plekhg6.png", Plekhg6.Egc_plot, width = 5, height = 5, dpi = 300)







######### H2-DMb1 ##########

# H2_DMb1.Tgc <- g.c.E.reversed[c("H2-DMb1","genotype")]
# colnames(H2_DMb1.Tgc)[1] <- "Expression_level"
# H2_DMb1.Tgc$genotype <- as.factor(H2_DMb1.Tgc$genotype)
# H2_DMb1.Tgc$Expression_level <- as.numeric(H2_DMb1.Tgc$Expression_level)
# png(file = "Plots/Liver_2wa_E_gc_H2-DMb1.png")
# H2_DMb1.Tgc_plot <- ggplot(H2_DMb1.Tgc,aes(x = genotype, y = Expression_level))+ 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of H2-DMb1 in Liver 2WA (gonad:chr) [E] ", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# H2_DMb1.Tgc_plot
# dev.off()

H2_DMb1.Egc <- g.c.E.reversed[,c("H2-DMb1","genotype")]
colnames(H2_DMb1.Egc)[1] <- "Expression_level"
H2_DMb1.Egc$genotype <- as.factor(H2_DMb1.Egc$genotype)
H2_DMb1.Egc$Expression_level <- as.numeric(H2_DMb1.Egc$Expression_level)

H2_DMb1.Egc_plot <- ggplot(H2_DMb1.Egc, aes(x= genotype, y= Expression_level)) + 
  geom_boxplot(outlier.shape = NA, fill = "#d04290") + 
  labs(title = "H2_DMb1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4, fill = "#d04290") 

## pairwise comparisons
H2_DMb1.Egc.stat <- H2_DMb1.Egc %>% 
  t_test(Expression_level ~ genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
H2_DMb1.Egc.stat

H2_DMb1.Egc.stat <- H2_DMb1.Egc.stat %>%
  add_xy_position(x = "genotype")

H2_DMb1.Egc_plot <- H2_DMb1.Egc_plot+
  stat_pvalue_manual(
    H2_DMb1.Egc.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Liver_2wa_E_gc/H2-DMb1.png", H2_DMb1.Egc_plot, width = 5, height = 5, dpi = 300)



###### arrange all plots together #######
liver.2wa.E.gc.plots <- ggarrange(H2_DMb1.Egc_plot, Plekhg6.Egc_plot,
                                  labels = c("A", "B"),
                                  ncol = 2, 
                                  nrow = 1)

ggsave("Plots_new/Liver_2wa_E_gc/Combined.png", liver.2wa.E.gc.plots, width = 10, height = 5, dpi = 600)
