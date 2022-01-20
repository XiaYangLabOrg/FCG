library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggeasy)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)
setwd("~/Desktop/Sex diff FCG mouse Jun2020-Jun2021/_submission/_Revision/Interactions_DEGs")

#####################===== T：gonad =====###################

Ing_3wa_T_GonadTreatment_DEGs <- read_excel("Data/Ing_3wa_T_GonadTreatment_DEGs.xlsx")
#(Ing_3wa_T_GonadTreatment_DEGs)[1] <- "Sample"
g.h.T <- Ing_3wa_T_GonadTreatment_DEGs %>% select(order(Ing_3wa_T_GonadTreatment_DEGs[3,]))
g.h.T.reversed <- t(g.h.T)
rownames(g.h.T.reversed)[1] <- "Sample"
colnames(g.h.T.reversed) <- g.h.T.reversed[1,]
#g.h.T.sample_name <- colnames(Ing_3wa_T_GonadTreatment_DEGs)[2:37]
g.h.T.reversed <- g.h.T.reversed[-1,]
g.h.T.samples <- rownames(g.h.T.reversed)[1:36]
g.h.T.samples <- as.data.frame(g.h.T.samples)
rownames(g.h.T.samples) <- g.h.T.samples[,1]
g.h.T.reversed <- g.h.T.reversed[-37,]
g.h.T.reversed <- cbind(g.h.T.samples,g.h.T.reversed)
g.h.T.reversed$genotype <- "temp"
g.h.T.reversed$genotype[1:7] <- "XXF"
g.h.T.reversed$genotype[8:17] <- "XXM"
g.h.T.reversed$genotype[18:27] <- "XYF"
g.h.T.reversed$genotype[28:36] <- "XYM"

colors_T <- c("B"="#ffffff" ,
              "T"="#2d4295")
############### Sephs1 #################

# Sephs1.Tgt <- g.t.T.reversed[,c("Sephs1","genotype","treatment")]
# colnames(Sephs1.Tgt)[1] <- "Expression_level"
# Sephs1.Tgt$genotype <- as.factor(Sephs1.Tgt$genotype)
# Sephs1.Tgt$Expression_level <- as.numeric(Sephs1.Tgt$Expression_level)
# png(file="Plots/Ing_3wa_T_gt_Sephs1.png")
# Sephs1.Tgt_plot <- ggplot(Sephs1.Tgt, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Sephs1 in Adipose (T+B:gonad)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Sephs1.Tgt_plot
# dev.off()

Sephs1.Tgh <- g.h.T.reversed[,c("Sephs1","genotype","treatment")]
colnames(Sephs1.Tgh)[1] <- "Expression_level"
Sephs1.Tgh$genotype <- as.factor(Sephs1.Tgh$genotype)
Sephs1.Tgh$Expression_level <- as.numeric(Sephs1.Tgh$Expression_level)

Sephs1.Tgh_plot <- ggplot(Sephs1.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Sephs1", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Sephs1.Tgh.stat <- Sephs1.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Sephs1.Tgh.stat

Sephs1.Tgh.stat <- Sephs1.Tgh.stat %>%
  add_xy_position(x = "genotype")

Sephs1.Tgh_plot <- Sephs1.Tgh_plot+
  stat_pvalue_manual(
    Sephs1.Tgh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_T_gh/Sephs1.png", Sephs1.Tgh_plot, width = 5, height = 5, dpi = 300)



######### Mpp5 #############

Mpp5.Tgh <- g.h.T.reversed[,c("Mpp5","genotype","treatment")]
colnames(Mpp5.Tgh)[1] <- "Expression_level"
Mpp5.Tgh$genotype <- as.factor(Mpp5.Tgh$genotype)
Mpp5.Tgh$Expression_level <- as.numeric(Mpp5.Tgh$Expression_level)

Mpp5.Tgh_plot <- ggplot(Mpp5.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Mpp5", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Mpp5.Tgh.stat <- Mpp5.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Mpp5.Tgh.stat

Mpp5.Tgh.stat <- Mpp5.Tgh.stat %>%
  add_xy_position(x = "genotype")

Mpp5.Tgh_plot <- Mpp5.Tgh_plot+
  stat_pvalue_manual(
    Mpp5.Tgh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_T_gh/Mpp5.png", Mpp5.Tgh_plot, width = 5, height = 5, dpi = 300)



#### arrange plots #####

ing.3wa.T.gt.plots <- ggarrange(Mpp5.Tgh_plot,Sephs1.Tgh_plot,
                                labels = c("A","B"),
                                ncol = 2, nrow = 1)

ggsave("Plots_new/Ing_3wa_T_gh/Combined.png", ing.3wa.T.gt.plots, width = 10, height = 5, dpi = 600)





#############======  T+B:chr  =====####################

# Mpp5
ggsave("Plots_new/Ing_3wa_T_ch/Mpp5.png", Mpp5.Tgh_plot, width = 5, height = 5, dpi = 300)


#############======  gonad:chr  =====####################

# Mpp5
ggsave("Plots_new/Ing_3wa_T_gc/Mpp5.png", Mpp5.Tgh_plot, width = 5, height = 5, dpi = 300)

###########======== chr:gonad:hormone =======#########

# Mpp5
ggsave("Plots_new/Ing_3wa_T_cgh/Mpp5.png", Mpp5.Tgh_plot, width = 5, height = 5, dpi = 300)



######## ====== E+B:gonad ======

Ing_3wa_E_GonadTreatment_DEGs <- read_excel("Data/Ing_3wa_E_GonadTreatment_DEGs.xlsx")
g.h.E <- Ing_3wa_E_GonadTreatment_DEGs %>% select(order(Ing_3wa_E_GonadTreatment_DEGs[32,]))
g.h.E.reversed <- t(g.h.E)
rownames(g.h.E.reversed)[1] <- "Sample"
colnames(g.h.E.reversed) <- g.h.E.reversed[1,]
g.h.E.reversed <- g.h.E.reversed[-1,]
g.h.E.samples <- rownames(g.h.E.reversed)
g.h.E.samples <- as.data.frame(g.h.E.samples)
rownames(g.h.E.samples) <- g.h.E.samples[,1]
g.h.E.reversed <- cbind(g.h.E.reversed,g.h.E.samples)
g.h.E.reversed$genotype <- "temp"
g.h.E.reversed$genotype[1:8] <- "XXF"
g.h.E.reversed$genotype[9:17] <- "XXM"
g.h.E.reversed$genotype[18:27] <- "XYF"
g.h.E.reversed$genotype[28:36] <- "XYM"

colors_E <- c("B"= "#ffffff",
              "E" = "#d04290")


######## Mpp5

# Mpp5.Egt <- g.h.E.reversed[,c("Mpp5","genotype","treatment")]
# colnames(Mpp5.Egt)[1] <- "Expression_level"
# Mpp5.Egt$genotype <- as.factor(Mpp5.Egt$genotype)
# Mpp5.Egt$Expression_level <- as.numeric(Mpp5.Egt$Expression_level)
# png(file="Plots/Ing_3wa_E_gt_Mpp5.png")
# Mpp5.Egt_plot <- ggplot(Mpp5.Egt, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Mpp5 in Adipose (E+B:gonad)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Mpp5.Egt_plot
# dev.off()

Mpp5.Egh <- g.h.E.reversed[,c("Mpp5","genotype","treatment")]
colnames(Mpp5.Egh)[1] <- "Expression_level"
Mpp5.Egh$genotype <- as.factor(Mpp5.Egh$genotype)
Mpp5.Egh$Expression_level <- as.numeric(Mpp5.Egh$Expression_level)

Mpp5.Egh_plot <- ggplot(Mpp5.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Mpp5", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Mpp5.Egh.stat <- Mpp5.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Mpp5.Egh.stat

Mpp5.Egh.stat <- Mpp5.Egh.stat %>%
  add_xy_position(x = "genotype")

Mpp5.Egh_plot <- Mpp5.Egh_plot+
  stat_pvalue_manual(
    Mpp5.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Mpp5.png", Mpp5.Egh_plot, width = 5, height = 5, dpi = 300)


# Acsl3
## The protein encoded by this gene is an isozyme of the long-chain fatty-acid-coenzyme A ligase family. Although differing in substrate specificity, subcellular localization, and tissue distribution, all isozymes of this family convert free long-chain fatty acids into fatty acyl-CoA esters, and thereby play a key role in lipid biosynthesis and fatty acid degradation. This isozyme is highly expressed in brain, and preferentially utilizes myristate, arachidonate, and eicosapentaenoate as substrates. The amino acid sequence of this isozyme is 92% identical to that of rat homolog. Two transcript variants encoding the same protein have been found for this gene. [provided by RefSeq, Jul 2008]


Acsl3.Egh <- g.h.E.reversed[,c("Acsl3","genotype","treatment")]
colnames(Acsl3.Egh)[1] <- "Expression_level"
Acsl3.Egh$genotype <- as.factor(Acsl3.Egh$genotype)
Acsl3.Egh$Expression_level <- as.numeric(Acsl3.Egh$Expression_level)

Acsl3.Egh_plot <- ggplot(Acsl3.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Acsl3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Acsl3.Egh.stat <- Acsl3.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Acsl3.Egh.stat

Acsl3.Egh.stat <- Acsl3.Egh.stat %>%
  add_xy_position(x = "genotype")

Acsl3.Egh_plot <- Acsl3.Egh_plot+
  stat_pvalue_manual(
    Acsl3.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Acsl3.png", Acsl3.Egh_plot, width = 5, height = 5, dpi = 300)

# Hsd11b1
## The protein encoded by this gene is a microsomal enzyme that catalyzes the conversion of the stress hormone cortisol to the inactive metabolite cortisone. In addition, the encoded protein can catalyze the reverse reaction, the conversion of cortisone to cortisol. Too much cortisol can lead to central obesity, and a particular variation in this gene has been associated with obesity and insulin resistance in children. Mutations in this gene and H6PD (hexose-6-phosphate dehydrogenase (glucose 1-dehydrogenase)) are the cause of cortisone reductase deficiency. Alternate splicing results in multiple transcript variants encoding the same protein.[provided by RefSeq, May 2011]
Hsd11b1.Egh <- g.h.E.reversed[,c("Hsd11b1","genotype","treatment")]
colnames(Hsd11b1.Egh)[1] <- "Expression_level"
Hsd11b1.Egh$genotype <- as.factor(Hsd11b1.Egh$genotype)
Hsd11b1.Egh$Expression_level <- as.numeric(Hsd11b1.Egh$Expression_level)

Hsd11b1.Egh_plot <- ggplot(Hsd11b1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Hsd11b1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Hsd11b1.Egh.stat <- Hsd11b1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Hsd11b1.Egh.stat

Hsd11b1.Egh.stat <- Hsd11b1.Egh.stat %>%
  add_xy_position(x = "genotype")

Hsd11b1.Egh_plot <- Hsd11b1.Egh_plot+
  stat_pvalue_manual(
    Hsd11b1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Hsd11b1.png", Hsd11b1.Egh_plot, width = 5, height = 5, dpi = 300)




###### Dnaic1

Dnaic1.Egh <- g.h.E.reversed[,c("Dnaic1","genotype","treatment")]
colnames(Dnaic1.Egh)[1] <- "Expression_level"
Dnaic1.Egh$genotype <- as.factor(Dnaic1.Egh$genotype)
Dnaic1.Egh$Expression_level <- as.numeric(Dnaic1.Egh$Expression_level)

Dnaic1.Egh_plot <- ggplot(Dnaic1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Dnaic1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Dnaic1.Egh.stat <- Dnaic1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Dnaic1.Egh.stat

Dnaic1.Egh.stat <- Dnaic1.Egh.stat %>%
  add_xy_position(x = "genotype")

Dnaic1.Egh_plot <- Dnaic1.Egh_plot+
  stat_pvalue_manual(
    Dnaic1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Dnaic1.png", Dnaic1.Egh_plot, width = 5, height = 5, dpi = 300)


##### Ptn
Ptn.Egh <- g.h.E.reversed[,c("Ptn","genotype","treatment")]
colnames(Ptn.Egh)[1] <- "Expression_level"
Ptn.Egh$genotype <- as.factor(Ptn.Egh$genotype)
Ptn.Egh$Expression_level <- as.numeric(Ptn.Egh$Expression_level)

Ptn.Egh_plot <- ggplot(Ptn.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Ptn", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Ptn.Egh.stat <- Ptn.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Ptn.Egh.stat

Ptn.Egh.stat <- Ptn.Egh.stat %>%
  add_xy_position(x = "genotype")

Ptn.Egh_plot <- Ptn.Egh_plot+
  stat_pvalue_manual(
    Ptn.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Ptn.png", Ptn.Egh_plot, width = 5, height = 5, dpi = 300)



##### Cited1
Cited1.Egh <- g.h.E.reversed[,c("Cited1","genotype","treatment")]
colnames(Cited1.Egh)[1] <- "Expression_level"
Cited1.Egh$genotype <- as.factor(Cited1.Egh$genotype)
Cited1.Egh$Expression_level <- as.numeric(Cited1.Egh$Expression_level)

Cited1.Egh_plot <- ggplot(Cited1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cited1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cited1.Egh.stat <- Cited1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cited1.Egh.stat

Cited1.Egh.stat <- Cited1.Egh.stat %>%
  add_xy_position(x = "genotype")

Cited1.Egh_plot <- Cited1.Egh_plot+
  stat_pvalue_manual(
    Cited1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Cited1.png", Cited1.Egh_plot, width = 5, height = 5, dpi = 300)



####### Ctns

Ctns.Egh <- g.h.E.reversed[,c("Ctns","genotype","treatment")]
colnames(Ctns.Egh)[1] <- "Expression_level"
Ctns.Egh$genotype <- as.factor(Ctns.Egh$genotype)
Ctns.Egh$Expression_level <- as.numeric(Ctns.Egh$Expression_level)

Ctns.Egh_plot <- ggplot(Ctns.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Ctns", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Ctns.Egh.stat <- Ctns.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Ctns.Egh.stat

Ctns.Egh.stat <- Ctns.Egh.stat %>%
  add_xy_position(x = "genotype")

Ctns.Egh_plot <- Ctns.Egh_plot+
  stat_pvalue_manual(
    Ctns.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Ctns.png", Ctns.Egh_plot, width = 5, height = 5, dpi = 300)




######### Slc2a3

Slc2a3.Egh <- g.h.E.reversed[,c("Slc2a3","genotype","treatment")]
colnames(Slc2a3.Egh)[1] <- "Expression_level"
Slc2a3.Egh$genotype <- as.factor(Slc2a3.Egh$genotype)
Slc2a3.Egh$Expression_level <- as.numeric(Slc2a3.Egh$Expression_level)

Slc2a3.Egh_plot <- ggplot(Slc2a3.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Slc2a3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Slc2a3.Egh.stat <- Slc2a3.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Slc2a3.Egh.stat

Slc2a3.Egh.stat <- Slc2a3.Egh.stat %>%
  add_xy_position(x = "genotype")

Slc2a3.Egh_plot <- Slc2a3.Egh_plot+
  stat_pvalue_manual(
    Slc2a3.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Slc2a3.png", Slc2a3.Egh_plot, width = 5, height = 5, dpi = 300)




#########  Tuft1

Tuft1.Egh <- g.h.E.reversed[,c("Tuft1","genotype","treatment")]
colnames(Tuft1.Egh)[1] <- "Expression_level"
Tuft1.Egh$genotype <- as.factor(Tuft1.Egh$genotype)
Tuft1.Egh$Expression_level <- as.numeric(Tuft1.Egh$Expression_level)

Tuft1.Egh_plot <- ggplot(Tuft1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Tuft1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Tuft1.Egh.stat <- Tuft1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Tuft1.Egh.stat

Tuft1.Egh.stat <- Tuft1.Egh.stat %>%
  add_xy_position(x = "genotype")

Tuft1.Egh_plot <- Tuft1.Egh_plot+
  stat_pvalue_manual(
    Tuft1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Tuft1.png", Tuft1.Egh_plot, width = 5, height = 5, dpi = 300)



###### Lamb3
Lamb3.Egh <- g.h.E.reversed[,c("Lamb3","genotype","treatment")]
colnames(Lamb3.Egh)[1] <- "Expression_level"
Lamb3.Egh$genotype <- as.factor(Lamb3.Egh$genotype)
Lamb3.Egh$Expression_level <- as.numeric(Lamb3.Egh$Expression_level)

Lamb3.Egh_plot <- ggplot(Lamb3.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Lamb3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Lamb3.Egh.stat <- Lamb3.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Lamb3.Egh.stat

Lamb3.Egh.stat <- Lamb3.Egh.stat %>%
  add_xy_position(x = "genotype")

Lamb3.Egh_plot <- Lamb3.Egh_plot+
  stat_pvalue_manual(
    Lamb3.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Lamb3.png", Lamb3.Egh_plot, width = 5, height = 5, dpi = 300)




###### Tmem86a

Tmem86a.Egh <- g.h.E.reversed[,c("Tmem86a","genotype","treatment")]
colnames(Tmem86a.Egh)[1] <- "Expression_level"
Tmem86a.Egh$genotype <- as.factor(Tmem86a.Egh$genotype)
Tmem86a.Egh$Expression_level <- as.numeric(Tmem86a.Egh$Expression_level)

Tmem86a.Egh_plot <- ggplot(Tmem86a.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Tmem86a", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Tmem86a.Egh.stat <- Tmem86a.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Tmem86a.Egh.stat

Tmem86a.Egh.stat <- Tmem86a.Egh.stat %>%
  add_xy_position(x = "genotype")

Tmem86a.Egh_plot <- Tmem86a.Egh_plot+
  stat_pvalue_manual(
    Tmem86a.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Tmem86a.png", Tmem86a.Egh_plot, width = 5, height = 5, dpi = 300)




######## Atp1b1
Atp1b1.Egh <- g.h.E.reversed[,c("Atp1b1","genotype","treatment")]
colnames(Atp1b1.Egh)[1] <- "Expression_level"
Atp1b1.Egh$genotype <- as.factor(Atp1b1.Egh$genotype)
Atp1b1.Egh$Expression_level <- as.numeric(Atp1b1.Egh$Expression_level)

Atp1b1.Egh_plot <- ggplot(Atp1b1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Atp1b1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Atp1b1.Egh.stat <- Atp1b1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Atp1b1.Egh.stat

Atp1b1.Egh.stat <- Atp1b1.Egh.stat %>%
  add_xy_position(x = "genotype")

Atp1b1.Egh_plot <- Atp1b1.Egh_plot+
  stat_pvalue_manual(
    Atp1b1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Atp1b1.png", Atp1b1.Egh_plot, width = 5, height = 5, dpi = 300)




######### Lmyc1

Lmyc1.Egh <- g.h.E.reversed[,c("Lmyc1","genotype","treatment")]
colnames(Lmyc1.Egh)[1] <- "Expression_level"
Lmyc1.Egh$genotype <- as.factor(Lmyc1.Egh$genotype)
Lmyc1.Egh$Expression_level <- as.numeric(Lmyc1.Egh$Expression_level)

Lmyc1.Egh_plot <- ggplot(Lmyc1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Lmyc1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Lmyc1.Egh.stat <- Lmyc1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Lmyc1.Egh.stat

Lmyc1.Egh.stat <- Lmyc1.Egh.stat %>%
  add_xy_position(x = "genotype")

Lmyc1.Egh_plot <- Lmyc1.Egh_plot+
  stat_pvalue_manual(
    Lmyc1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Lmyc1.png", Lmyc1.Egh_plot, width = 5, height = 5, dpi = 300)




######## Mfsd2a

Mfsd2a.Egh <- g.h.E.reversed[,c("Mfsd2a","genotype","treatment")]
colnames(Mfsd2a.Egh)[1] <- "Expression_level"
Mfsd2a.Egh$genotype <- as.factor(Mfsd2a.Egh$genotype)
Mfsd2a.Egh$Expression_level <- as.numeric(Mfsd2a.Egh$Expression_level)

Mfsd2a.Egh_plot <- ggplot(Mfsd2a.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Mfsd2a", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Mfsd2a.Egh.stat <- Mfsd2a.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Mfsd2a.Egh.stat

Mfsd2a.Egh.stat <- Mfsd2a.Egh.stat %>%
  add_xy_position(x = "genotype")

Mfsd2a.Egh_plot <- Mfsd2a.Egh_plot+
  stat_pvalue_manual(
    Mfsd2a.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Mfsd2a.png", Mfsd2a.Egh_plot, width = 5, height = 5, dpi = 300)


######### Sytl2
Sytl2.Egh <- g.h.E.reversed[,c("Sytl2","genotype","treatment")]
colnames(Sytl2.Egh)[1] <- "Expression_level"
Sytl2.Egh$genotype <- as.factor(Sytl2.Egh$genotype)
Sytl2.Egh$Expression_level <- as.numeric(Sytl2.Egh$Expression_level)

Sytl2.Egh_plot <- ggplot(Sytl2.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Sytl2", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Sytl2.Egh.stat <- Sytl2.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Sytl2.Egh.stat

Sytl2.Egh.stat <- Sytl2.Egh.stat %>%
  add_xy_position(x = "genotype")

Sytl2.Egh_plot <- Sytl2.Egh_plot+
  stat_pvalue_manual(
    Sytl2.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Sytl2.png", Sytl2.Egh_plot, width = 5, height = 5, dpi = 300)




######## Tmed3
Tmed3.Egh <- g.h.E.reversed[,c("Tmed3","genotype","treatment")]
colnames(Tmed3.Egh)[1] <- "Expression_level"
Tmed3.Egh$genotype <- as.factor(Tmed3.Egh$genotype)
Tmed3.Egh$Expression_level <- as.numeric(Tmed3.Egh$Expression_level)

Tmed3.Egh_plot <- ggplot(Tmed3.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Tmed3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Tmed3.Egh.stat <- Tmed3.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Tmed3.Egh.stat

Tmed3.Egh.stat <- Tmed3.Egh.stat %>%
  add_xy_position(x = "genotype")

Tmed3.Egh_plot <- Tmed3.Egh_plot+
  stat_pvalue_manual(
    Tmed3.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Tmed3.png", Tmed3.Egh_plot, width = 5, height = 5, dpi = 300)




######## Nos3as
Nos3as.Egh <- g.h.E.reversed[,c("Nos3as","genotype","treatment")]
colnames(Nos3as.Egh)[1] <- "Expression_level"
Nos3as.Egh$genotype <- as.factor(Nos3as.Egh$genotype)
Nos3as.Egh$Expression_level <- as.numeric(Nos3as.Egh$Expression_level)

Nos3as.Egh_plot <- ggplot(Nos3as.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Nos3as", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Nos3as.Egh.stat <- Nos3as.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Nos3as.Egh.stat

Nos3as.Egh.stat <- Nos3as.Egh.stat %>%
  add_xy_position(x = "genotype")

Nos3as.Egh_plot <- Nos3as.Egh_plot+
  stat_pvalue_manual(
    Nos3as.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Nos3as.png", Nos3as.Egh_plot, width = 5, height = 5, dpi = 300)




####### Pxylp1
Pxylp1.Egh <- g.h.E.reversed[,c("Pxylp1","genotype","treatment")]
colnames(Pxylp1.Egh)[1] <- "Expression_level"
Pxylp1.Egh$genotype <- as.factor(Pxylp1.Egh$genotype)
Pxylp1.Egh$Expression_level <- as.numeric(Pxylp1.Egh$Expression_level)

Pxylp1.Egh_plot <- ggplot(Pxylp1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Pxylp1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Pxylp1.Egh.stat <- Pxylp1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Pxylp1.Egh.stat

Pxylp1.Egh.stat <- Pxylp1.Egh.stat %>%
  add_xy_position(x = "genotype")

Pxylp1.Egh_plot <- Pxylp1.Egh_plot+
  stat_pvalue_manual(
    Pxylp1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Pxylp1.png", Pxylp1.Egh_plot, width = 5, height = 5, dpi = 300)




####### Adra2a
Adra2a.Egh <- g.h.E.reversed[,c("Adra2a","genotype","treatment")]
colnames(Adra2a.Egh)[1] <- "Expression_level"
Adra2a.Egh$genotype <- as.factor(Adra2a.Egh$genotype)
Adra2a.Egh$Expression_level <- as.numeric(Adra2a.Egh$Expression_level)

Adra2a.Egh_plot <- ggplot(Adra2a.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Adra2a", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Adra2a.Egh.stat <- Adra2a.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Adra2a.Egh.stat

Adra2a.Egh.stat <- Adra2a.Egh.stat %>%
  add_xy_position(x = "genotype")

Adra2a.Egh_plot <- Adra2a.Egh_plot+
  stat_pvalue_manual(
    Adra2a.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Adra2a.png", Adra2a.Egh_plot, width = 5, height = 5, dpi = 300)



###### Trp63
Trp63.Egh <- g.h.E.reversed[,c("Trp63","genotype","treatment")]
colnames(Trp63.Egh)[1] <- "Expression_level"
Trp63.Egh$genotype <- as.factor(Trp63.Egh$genotype)
Trp63.Egh$Expression_level <- as.numeric(Trp63.Egh$Expression_level)

Trp63.Egh_plot <- ggplot(Trp63.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Trp63", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Trp63.Egh.stat <- Trp63.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Trp63.Egh.stat

Trp63.Egh.stat <- Trp63.Egh.stat %>%
  add_xy_position(x = "genotype")

Trp63.Egh_plot <- Trp63.Egh_plot+
  stat_pvalue_manual(
    Trp63.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Trp63.png", Trp63.Egh_plot, width = 5, height = 5, dpi = 300)




########### Dnajc12
Dnajc12.Egh <- g.h.E.reversed[,c("Dnajc12","genotype","treatment")]
colnames(Dnajc12.Egh)[1] <- "Expression_level"
Dnajc12.Egh$genotype <- as.factor(Dnajc12.Egh$genotype)
Dnajc12.Egh$Expression_level <- as.numeric(Dnajc12.Egh$Expression_level)

Dnajc12.Egh_plot <- ggplot(Dnajc12.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Dnajc12", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Dnajc12.Egh.stat <- Dnajc12.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Dnajc12.Egh.stat

Dnajc12.Egh.stat <- Dnajc12.Egh.stat %>%
  add_xy_position(x = "genotype")

Dnajc12.Egh_plot <- Dnajc12.Egh_plot+
  stat_pvalue_manual(
    Dnajc12.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Dnajc12.png", Dnajc12.Egh_plot, width = 5, height = 5, dpi = 300)




######### Gas6
Gas6.Egh <- g.h.E.reversed[,c("Gas6","genotype","treatment")]
colnames(Gas6.Egh)[1] <- "Expression_level"
Gas6.Egh$genotype <- as.factor(Gas6.Egh$genotype)
Gas6.Egh$Expression_level <- as.numeric(Gas6.Egh$Expression_level)

Gas6.Egh_plot <- ggplot(Gas6.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Gas6", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Gas6.Egh.stat <- Gas6.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Gas6.Egh.stat

Gas6.Egh.stat <- Gas6.Egh.stat %>%
  add_xy_position(x = "genotype")

Gas6.Egh_plot <- Gas6.Egh_plot+
  stat_pvalue_manual(
    Gas6.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Gas6.png", Gas6.Egh_plot, width = 5, height = 5, dpi = 300)



########## Anxa8
Anxa8.Egh <- g.h.E.reversed[,c("Anxa8","genotype","treatment")]
colnames(Anxa8.Egh)[1] <- "Expression_level"
Anxa8.Egh$genotype <- as.factor(Anxa8.Egh$genotype)
Anxa8.Egh$Expression_level <- as.numeric(Anxa8.Egh$Expression_level)

Anxa8.Egh_plot <- ggplot(Anxa8.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Anxa8", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Anxa8.Egh.stat <- Anxa8.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Anxa8.Egh.stat

Anxa8.Egh.stat <- Anxa8.Egh.stat %>%
  add_xy_position(x = "genotype")

Anxa8.Egh_plot <- Anxa8.Egh_plot+
  stat_pvalue_manual(
    Anxa8.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Anxa8.png", Anxa8.Egh_plot, width = 5, height = 5, dpi = 300)



####### Myl12a

Myl12a.Egh <- g.h.E.reversed[,c("Myl12a","genotype","treatment")]
colnames(Myl12a.Egh)[1] <- "Expression_level"
Myl12a.Egh$genotype <- as.factor(Myl12a.Egh$genotype)
Myl12a.Egh$Expression_level <- as.numeric(Myl12a.Egh$Expression_level)

Myl12a.Egh_plot <- ggplot(Myl12a.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Myl12a", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Myl12a.Egh.stat <- Myl12a.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Myl12a.Egh.stat

Myl12a.Egh.stat <- Myl12a.Egh.stat %>%
  add_xy_position(x = "genotype")

Myl12a.Egh_plot <- Myl12a.Egh_plot+
  stat_pvalue_manual(
    Myl12a.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Myl12a.png", Myl12a.Egh_plot, width = 5, height = 5, dpi = 300)




######## Lrrc26
Lrrc26.Egh <- g.h.E.reversed[,c("Lrrc26","genotype","treatment")]
colnames(Lrrc26.Egh)[1] <- "Expression_level"
Lrrc26.Egh$genotype <- as.factor(Lrrc26.Egh$genotype)
Lrrc26.Egh$Expression_level <- as.numeric(Lrrc26.Egh$Expression_level)

Lrrc26.Egh_plot <- ggplot(Lrrc26.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Lrrc26", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Lrrc26.Egh.stat <- Lrrc26.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Lrrc26.Egh.stat

Lrrc26.Egh.stat <- Lrrc26.Egh.stat %>%
  add_xy_position(x = "genotype")

Lrrc26.Egh_plot <- Lrrc26.Egh_plot+
  stat_pvalue_manual(
    Lrrc26.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Lrrc26.png", Lrrc26.Egh_plot, width = 5, height = 5, dpi = 300)




######## Tapt1

Tapt1.Egh <- g.h.E.reversed[,c("Tapt1","genotype","treatment")]
colnames(Tapt1.Egh)[1] <- "Expression_level"
Tapt1.Egh$genotype <- as.factor(Tapt1.Egh$genotype)
Tapt1.Egh$Expression_level <- as.numeric(Tapt1.Egh$Expression_level)

Tapt1.Egh_plot <- ggplot(Tapt1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Tapt1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Tapt1.Egh.stat <- Tapt1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Tapt1.Egh.stat

Tapt1.Egh.stat <- Tapt1.Egh.stat %>%
  add_xy_position(x = "genotype")

Tapt1.Egh_plot <- Tapt1.Egh_plot+
  stat_pvalue_manual(
    Tapt1.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Tapt1.png", Tapt1.Egh_plot, width = 5, height = 5, dpi = 300)



########### Anln
Anln.Egh <- g.h.E.reversed[,c("Anln","genotype","treatment")]
colnames(Anln.Egh)[1] <- "Expression_level"
Anln.Egh$genotype <- as.factor(Anln.Egh$genotype)
Anln.Egh$Expression_level <- as.numeric(Anln.Egh$Expression_level)

Anln.Egh_plot <- ggplot(Anln.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Anln", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Anln.Egh.stat <- Anln.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Anln.Egh.stat

Anln.Egh.stat <- Anln.Egh.stat %>%
  add_xy_position(x = "genotype")

Anln.Egh_plot <- Anln.Egh_plot+
  stat_pvalue_manual(
    Anln.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Anln.png", Anln.Egh_plot, width = 5, height = 5, dpi = 300)



######### Stc2

Stc2.Egh <- g.h.E.reversed[,c("Stc2","genotype","treatment")]
colnames(Stc2.Egh)[1] <- "Expression_level"
Stc2.Egh$genotype <- as.factor(Stc2.Egh$genotype)
Stc2.Egh$Expression_level <- as.numeric(Stc2.Egh$Expression_level)

Stc2.Egh_plot <- ggplot(Stc2.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Stc2", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Stc2.Egh.stat <- Stc2.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Stc2.Egh.stat

Stc2.Egh.stat <- Stc2.Egh.stat %>%
  add_xy_position(x = "genotype")

Stc2.Egh_plot <- Stc2.Egh_plot+
  stat_pvalue_manual(
    Stc2.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Stc2.png", Stc2.Egh_plot, width = 5, height = 5, dpi = 300)



###### Tlcd3a
Tlcd3a.Egh <- g.h.E.reversed[,c("Tlcd3a","genotype","treatment")]
colnames(Tlcd3a.Egh)[1] <- "Expression_level"
Tlcd3a.Egh$genotype <- as.factor(Tlcd3a.Egh$genotype)
Tlcd3a.Egh$Expression_level <- as.numeric(Tlcd3a.Egh$Expression_level)

Tlcd3a.Egh_plot <- ggplot(Tlcd3a.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Tlcd3a", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Tlcd3a.Egh.stat <- Tlcd3a.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Tlcd3a.Egh.stat

Tlcd3a.Egh.stat <- Tlcd3a.Egh.stat %>%
  add_xy_position(x = "genotype")

Tlcd3a.Egh_plot <- Tlcd3a.Egh_plot+
  stat_pvalue_manual(
    Tlcd3a.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Tlcd3a.png", Tlcd3a.Egh_plot, width = 5, height = 5, dpi = 300)




###### S100a14

S100a14.Egh <- g.h.E.reversed[,c("S100a14","genotype","treatment")]
colnames(S100a14.Egh)[1] <- "Expression_level"
S100a14.Egh$genotype <- as.factor(S100a14.Egh$genotype)
S100a14.Egh$Expression_level <- as.numeric(S100a14.Egh$Expression_level)

S100a14.Egh_plot <- ggplot(S100a14.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "S100a14", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
S100a14.Egh.stat <- S100a14.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
S100a14.Egh.stat

S100a14.Egh.stat <- S100a14.Egh.stat %>%
  add_xy_position(x = "genotype")

S100a14.Egh_plot <- S100a14.Egh_plot+
  stat_pvalue_manual(
    S100a14.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/S100a14.png", S100a14.Egh_plot, width = 5, height = 5, dpi = 300)




##### Ier3

Ier3.Egh <- g.h.E.reversed[,c("Ier3","genotype","treatment")]
colnames(Ier3.Egh)[1] <- "Expression_level"
Ier3.Egh$genotype <- as.factor(Ier3.Egh$genotype)
Ier3.Egh$Expression_level <- as.numeric(Ier3.Egh$Expression_level)

Ier3.Egh_plot <- ggplot(Ier3.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Ier3", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Ier3.Egh.stat <- Ier3.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Ier3.Egh.stat

Ier3.Egh.stat <- Ier3.Egh.stat %>%
  add_xy_position(x = "genotype")

Ier3.Egh_plot <- Ier3.Egh_plot+
  stat_pvalue_manual(
    Ier3.Egh.stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Ing_3wa_E_gh/Ier3.png", Ier3.Egh_plot, width = 5, height = 5, dpi = 300)


##### arrange

Ing.3wa.E.gh.plots <- ggarrange(Dnaic1.Egh_plot, Ptn.Egh_plot, Cited1.Egh_plot, Ctns.Egh_plot, Slc2a3.Egh_plot, Tuft1.Egh_plot, Lamb3.Egh_plot, Tmem86a.Egh_plot, Atp1b1.Egh_plot, Lmyc1.Egh_plot, Mfsd2a.Egh_plot, Sytl2.Egh_plot, Tmed3.Egh_plot, Nos3as.Egh_plot, Pxylp1.Egh_plot, Adra2a.Egh_plot, Trp63.Egh_plot, Hsd11b1.Egh_plot, Dnajc12.Egh_plot, Gas6.Egh_plot, Anxa8.Egh_plot, Myl12a.Egh_plot, Mpp5.Egh_plot, Lrrc26.Egh_plot, Tapt1.Egh_plot, Anln.Egh_plot, Stc2.Egh_plot, Tlcd3a.Egh_plot, Acsl3.Egh_plot, S100a14.Egh_plot, Ier3.Egh_plot,
                               ncol = 11, nrow = 3)

ggsave("Plots_new/Ing_3wa_E_gh/Combined.png",Ing.3wa.E.gh.plots,width = 50, height = 15, limitsize = FALSE, dpi = 600)






############## chr: E+B 

ggsave("Plots_new/Ing_3wa_E_gc/Dnaic1.png", Dnaic1.Egh_plot, width = 5, height = 5, dpi = 300)



