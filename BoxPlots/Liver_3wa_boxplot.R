library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggeasy)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)


##################========= T:gonad ==============#####################
colors_T <- c("B"="#ffffff" ,
              "T"="#2d4295")

liver_3wa_T_GonadTreatment_DEGs <- read_excel("Data/liver_3wa_T_GonadTreatment_DEGs.xlsx")
colnames(liver_3wa_T_GonadTreatment_DEGs)[1] <- "Sample"
g.t <- liver_3wa_T_GonadTreatment_DEGs %>% select(order(liver_3wa_T_GonadTreatment_DEGs[11,]))
g.t <- g.t[,-2]
g.h.reversed <- t(g.t)
colnames(g.h.reversed) <- g.h.reversed[1,]
sample_name <- colnames(liver_3wa_T_GonadTreatment_DEGs)[2:41]
g.h.reversed <- g.h.reversed[-1,]
samples <- rownames(g.h.reversed)
samples <- as.data.frame(samples)
rownames(samples) <- samples[,1]
g.h.reversed <- cbind(samples, g.h.reversed)
g.h.reversed$genotype <- "temp"
g.h.reversed$genotype[1:10] <- "XXF"
g.h.reversed$genotype[11:20] <- "XXM"
g.h.reversed$genotype[21:30] <- "XYF"
g.h.reversed$genotype[31:40] <- "XYM"


####### Sult3a1 ########
Sult3a1.Tgh <- g.h.reversed[,c("Sult3a1","genotype","treatment")]
colnames(Sult3a1.Tgh)[1] <- "Expression_level"
Sult3a1.Tgh$genotype <- as.factor(Sult3a1.Tgh$genotype)
Sult3a1.Tgh$Expression_level <- as.numeric(Sult3a1.Tgh$Expression_level)

Sult3a1.Tgh_plot <- ggplot(Sult3a1.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Sult3a1", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Sult3a1.Tgh.stat <- Sult3a1.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Sult3a1.Tgh.stat

Sult3a1.Tgh.stat <- Sult3a1.Tgh.stat %>%
  add_xy_position(x = "genotype")

Sult3a1.Tgh_plot <- Sult3a1.Tgh_plot+
  stat_pvalue_manual(
    Sult3a1.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Sult3a1.png", Sult3a1.Tgh_plot, width = 5, height = 5, dpi = 300)





######## BC027556 ##########

## is also known as Lcn13 : Predicted to enable small molecule binding activity. Acts upstream of or within several processes, including insulin receptor signaling pathway; negative regulation of gluconeogenesis; and regulation of lipid metabolic process. Located in extracellular space.Orthologous to human OBP2A (odorant binding protein 2A. 

# Lcn13.Tgh <- g.h.reversed[,c("BC027556","genotype","treatment")]
# colnames(Lcn13.Tgh)[1] <- "Expression_level"
# Lcn13.Tgh$genotype <- as.factor(Lcn13.Tgh$genotype)
# Lcn13.Tgh$Expression_level <- as.numeric(Lcn13.Tgh$Expression_level)
# png(file="Plots/liver_3wa_T_gh_Lcn13.png")
# Lcn13.Tgh_plot <- ggplot(Lcn13.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Lcn13 in Liver (T+B:Gonad)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Lcn13.Tgh_plot
# dev.off()


Lcn13.Tgh <- g.h.reversed[,c("BC027556","genotype","treatment")]
colnames(Lcn13.Tgh)[1] <- "Expression_level"
Lcn13.Tgh$genotype <- as.factor(Lcn13.Tgh$genotype)
Lcn13.Tgh$Expression_level <- as.numeric(Lcn13.Tgh$Expression_level)

Lcn13.Tgh_plot <- ggplot(Lcn13.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Lcn13", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Lcn13.Tgh.stat <- Lcn13.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Lcn13.Tgh.stat

Lcn13.Tgh.stat <- Lcn13.Tgh.stat %>%
  add_xy_position(x = "genotype")

Lcn13.Tgh_plot <- Lcn13.Tgh_plot+
  stat_pvalue_manual(
    Lcn13.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Lcn13.png", Lcn13.Tgh_plot, width = 5, height = 5, dpi = 300)




########### Cyp17a1########### 

# This gene encodes a member of the cytochrome P450 superfamily of enzymes. The cytochrome P450 proteins are monooxygenases which catalyze many reactions involved in drug metabolism and synthesis of cholesterol, steroids and other lipids. This protein localizes to the endoplasmic reticulum. It has both 17alpha-hydroxylase and 17,20-lyase activities and is a key enzyme in the steroidogenic pathway that produces progestins, mineralocorticoids, glucocorticoids, androgens, and estrogens.

# Cyp17a1.Tgh <- g.h.reversed[,c("Cyp17a1","group","treatment","Gonad")]
# colnames(Cyp17a1.Tgh)[1] <- "Expression_level"
# Cyp17a1.Tgh$group <- as.factor(Cyp17a1.Tgh$group)
# Cyp17a1.Tgh$Expression_level <- as.numeric(Cyp17a1.Tgh$Expression_level)
# Cyp17a1.Tgh_plot <- ggplot(Cyp17a1.Tgh, aes(x= Gonad, y= Expression_level, fill = group)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Cyp17a1 in Liver (T)", y = "Expression Level", x = "Gonadal Sex")+
#   scale_fill_manual(values= colors,limits= c("XX_B","XX_T","XY-_B","XY-_T","XXSry_B","XXSry_T","XY-Sry_B","XY-Sry_T"))+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Cyp17a1.Tgh_plot

# Cyp17a1.Tgh.Tgh <- g.h.reversed[,c("Cyp17a1","genotype","treatment")]
# colnames(Cyp17a1.Tgh)[1] <- "Expression_level"
# Cyp17a1.Tgh$genotype <- as.factor(Cyp17a1.Tgh$genotype)
# Cyp17a1.Tgh$Expression_level <- as.numeric(Cyp17a1.Tgh$Expression_level)
# png(file="Plots/liver_3wa_T_gh_Cyp17a1.png")
# Cyp17a1.Tgh_plot <- ggplot(Cyp17a1.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Cyp17a1 in Liver (T+B:Goand)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Cyp17a1.Tgh_plot
# dev.off()

Cyp17a1.Tgh <- g.h.reversed[,c("Cyp17a1","genotype","treatment")]
colnames(Cyp17a1.Tgh)[1] <- "Expression_level"
Cyp17a1.Tgh$genotype <- as.factor(Cyp17a1.Tgh$genotype)
Cyp17a1.Tgh$Expression_level <- as.numeric(Cyp17a1.Tgh$Expression_level)

Cyp17a1.Tgh_plot <- ggplot(Cyp17a1.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cyp17a1", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cyp17a1.Tgh.stat <- Cyp17a1.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cyp17a1.Tgh.stat

Cyp17a1.Tgh.stat <- Cyp17a1.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cyp17a1.Tgh_plot <- Cyp17a1.Tgh_plot+
  stat_pvalue_manual(
    Cyp17a1.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cyp17a1.png", Cyp17a1.Tgh_plot, width = 5, height = 5, dpi = 300)



########### Cyp3a41 ##############
Cyp3a41.Tgh <- g.h.reversed[,c("Cyp3a41","genotype","treatment")]
colnames(Cyp3a41.Tgh)[1] <- "Expression_level"
Cyp3a41.Tgh$genotype <- as.factor(Cyp3a41.Tgh$genotype)
Cyp3a41.Tgh$Expression_level <- as.numeric(Cyp3a41.Tgh$Expression_level)

Cyp3a41.Tgh_plot <- ggplot(Cyp3a41.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cyp3a41", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cyp3a41.Tgh.stat <- Cyp3a41.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cyp3a41.Tgh.stat

Cyp3a41.Tgh.stat <- Cyp3a41.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cyp3a41.Tgh_plot <- Cyp3a41.Tgh_plot+
  stat_pvalue_manual(
    Cyp3a41.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cyp3a41.png", Cyp3a41.Tgh_plot, width = 5, height = 5, dpi = 300)




####### Cyp2d9 ########
Cyp2d9.Tgh <- g.h.reversed[,c("Cyp2d9","genotype","treatment")]
colnames(Cyp2d9.Tgh)[1] <- "Expression_level"
Cyp2d9.Tgh$genotype <- as.factor(Cyp2d9.Tgh$genotype)
Cyp2d9.Tgh$Expression_level <- as.numeric(Cyp2d9.Tgh$Expression_level)

Cyp2d9.Tgh_plot <- ggplot(Cyp2d9.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cyp2d9", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"), 
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cyp2d9.Tgh.stat <- Cyp2d9.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cyp2d9.Tgh.stat

Cyp2d9.Tgh.stat <- Cyp2d9.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cyp2d9.Tgh_plot <- Cyp2d9.Tgh_plot+
  stat_pvalue_manual(
    Cyp2d9.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cyp2d9.png", Cyp2d9.Tgh_plot, width = 5, height = 5, dpi = 300)






######### Igfbp2 ##########

Igfbp2.Tgh <- g.h.reversed[,c("Igfbp2","genotype","treatment")]
colnames(Igfbp2.Tgh)[1] <- "Expression_level"
Igfbp2.Tgh$genotype <- as.factor(Igfbp2.Tgh$genotype)
Igfbp2.Tgh$Expression_level <- as.numeric(Igfbp2.Tgh$Expression_level)

Igfbp2.Tgh_plot <- ggplot(Igfbp2.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Igfbp2", y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Igfbp2.Tgh.stat <- Igfbp2.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Igfbp2.Tgh.stat

Igfbp2.Tgh.stat <- Igfbp2.Tgh.stat %>%
  add_xy_position(x = "genotype")

Igfbp2.Tgh_plot <- Igfbp2.Tgh_plot+
  stat_pvalue_manual(
    Igfbp2.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Igfbp2.png", Igfbp2.Tgh_plot, width = 5, height = 5, dpi = 300)




###### Cml4 #########

Cml4.Tgh <- g.h.reversed[,c("Cml4","genotype","treatment")]
colnames(Cml4.Tgh)[1] <- "Expression_level"
Cml4.Tgh$genotype <- as.factor(Cml4.Tgh$genotype)
Cml4.Tgh$Expression_level <- as.numeric(Cml4.Tgh$Expression_level)

Cml4.Tgh_plot <- ggplot(Cml4.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cml4", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cml4.Tgh.stat <- Cml4.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cml4.Tgh.stat

Cml4.Tgh.stat <- Cml4.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cml4.Tgh_plot <- Cml4.Tgh_plot+
  stat_pvalue_manual(
    Cml4.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cml4.png", Cml4.Tgh_plot, width = 5, height = 5, dpi = 300)





######## Cyp4a12 ###########

Cyp4a12.Tgh <- g.h.reversed[,c("Cyp4a12","genotype","treatment")]
colnames(Cyp4a12.Tgh)[1] <- "Expression_level"
Cyp4a12.Tgh$genotype <- as.factor(Cyp4a12.Tgh$genotype)
Cyp4a12.Tgh$Expression_level <- as.numeric(Cyp4a12.Tgh$Expression_level)

Cyp4a12.Tgh_plot <- ggplot(Cyp4a12.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cyp4a12", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cyp4a12.Tgh.stat <- Cyp4a12.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cyp4a12.Tgh.stat

Cyp4a12.Tgh.stat <- Cyp4a12.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cyp4a12.Tgh_plot <- Cyp4a12.Tgh_plot+
  stat_pvalue_manual(
    Cyp4a12.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cyp4a12.png", Cyp4a12.Tgh_plot, width = 5, height = 5, dpi = 300)



######## Susd4 (E430021N18Rik) #############

Susd4.Tgh <- g.h.reversed[,c("E430021N18Rik","genotype","treatment")]
colnames(Susd4.Tgh)[1] <- "Expression_level"
Susd4.Tgh$genotype <- as.factor(Susd4.Tgh$genotype)
Susd4.Tgh$Expression_level <- as.numeric(Susd4.Tgh$Expression_level)

Susd4.Tgh_plot <- ggplot(Susd4.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Susd4", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Susd4.Tgh.stat <- Susd4.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Susd4.Tgh.stat

Susd4.Tgh.stat <- Susd4.Tgh.stat %>%
  add_xy_position(x = "genotype")

Susd4.Tgh_plot <- Susd4.Tgh_plot+
  stat_pvalue_manual(
    Susd4.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Susd4.png", Susd4.Tgh_plot, width = 5, height = 5, dpi = 300)





######### Cxcl9 #########

Cxcl9.Tgh <- g.h.reversed[,c("Cxcl9","genotype","treatment")]
colnames(Cxcl9.Tgh)[1] <- "Expression_level"
Cxcl9.Tgh$genotype <- as.factor(Cxcl9.Tgh$genotype)
Cxcl9.Tgh$Expression_level <- as.numeric(Cxcl9.Tgh$Expression_level)

Cxcl9.Tgh_plot <- ggplot(Cxcl9.Tgh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Cxcl9", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_T,limits= c("T","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Cxcl9.Tgh.stat <- Cxcl9.Tgh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Cxcl9.Tgh.stat

Cxcl9.Tgh.stat <- Cxcl9.Tgh.stat %>%
  add_xy_position(x = "genotype")

Cxcl9.Tgh_plot <- Cxcl9.Tgh_plot+
  stat_pvalue_manual(
    Cxcl9.Tgh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/liver_3wa_T_gh_Cxcl9.png", Cxcl9.Tgh_plot, width = 5, height = 5, dpi = 300)





##### arrange all DEGs into one figure ##########

liver.3wa.T.gh.plots <- ggarrange(Cyp3a41.Tgh_plot, Sult3a1.Tgh_plot, Cyp2d9.Tgh_plot, Lcn13.Tgh_plot, Igfbp2.Tgh_plot, Cyp17a1.Tgh_plot, Cml4.Tgh_plot, Cyp4a12.Tgh_plot, Susd4.Tgh_plot, Cxcl9.Tgh_plot, 
                                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                                  ncol = 2, 
                                  nrow = 5)

ggsave("Plots_new/Liver_3wa_T_gh/Combined.png", liver.3wa.T.gh.plots, width = 15, height = 20, dpi = 600)

















#########################========= E+B:gonad ===========################################

colors_E <- c("B"= "#ffffff",
              "E" = "#d04290")

liver_3wa_E_GonadTreatment_DEGs <- read_excel("Data/Liver_3wa_E_GonadTreatment_DEGs.xlsx")
colnames(liver_3wa_E_GonadTreatment_DEGs)[1] <- "Sample"
g.h.E <- liver_3wa_E_GonadTreatment_DEGs %>% select(order(liver_3wa_E_GonadTreatment_DEGs[3,]))
g.h.E.reversed <- t(g.h.E)
rownames(g.h.E.reversed)[1] <- "Sample"
colnames(g.h.E.reversed) <- g.h.E.reversed[1,]
g.h.E.sample_name <- colnames(liver_3wa_E_GonadTreatment_DEGs)[2:41]
g.h.E.reversed <- g.h.E.reversed[-1,]
g.h.E.samples <- rownames(g.h.E.reversed)
g.h.E.samples <- as.data.frame(g.h.E.samples)
rownames(g.h.E.samples) <- g.h.E.samples[,1]
g.h.E.reversed <- cbind(g.h.E.samples, g.h.E.reversed)
g.h.E.reversed$genotype <- "temp"
g.h.E.reversed$genotype[1:10] <- "XXF"
g.h.E.reversed$genotype[11:20] <- "XXM"
g.h.E.reversed$genotype[21:30] <- "XYF"
g.h.E.reversed$genotype[31:40] <- "XYM"


###### Sult3a1 #######

# Sult3a1.Egh <- g.h.E.reversed[,c("Sult3a1","genotype","treatment")]
# colnames(Sult3a1.Egh)[1] <- "Expression_level"
# Sult3a1.Egh$genotype <- as.factor(Sult3a1.Egh$genotype)
# Sult3a1.Egh$Expression_level <- as.numeric(Sult3a1.Egh$Expression_level)
# png(file="Plots/liver_3wa_E_gh_Sult3a1.png")
# Sult3a1.Egh_plot <- ggplot(Sult3a1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Sult3a1 in Liver (E+B:gonad)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Sult3a1.Egh_plot
# dev.off()


Sult3a1.Egh <- g.h.E.reversed[,c("Sult3a1","genotype","treatment")]
colnames(Sult3a1.Egh)[1] <- "Expression_level"
Sult3a1.Egh$genotype <- as.factor(Sult3a1.Egh$genotype)
Sult3a1.Egh$Expression_level <- as.numeric(Sult3a1.Egh$Expression_level)

Sult3a1.Egh_plot <- ggplot(Sult3a1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Sult3a1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Sult3a1.Egh.stat <- Sult3a1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Sult3a1.Egh.stat

Sult3a1.Egh.stat <- Sult3a1.Egh.stat %>%
  add_xy_position(x = "genotype")

Sult3a1.Egh_plot <- Sult3a1.Egh_plot+
  stat_pvalue_manual(
    Sult3a1.Egh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Liver_3wa_E_gh/Sult3a1.png", Sult3a1.Egh_plot, width = 5, height = 5, dpi = 300)

######## Slc11a1 ##########

# Slc11a1.Egh <- g.h.E.reversed[,c("Slc11a1","genotype","treatment")]
# colnames(Slc11a1.Egh)[1] <- "Expression_level"
# Slc11a1.Egh$genotype <- as.factor(Slc11a1.Egh$genotype)
# Slc11a1.Egh$Expression_level <- as.numeric(Slc11a1.Egh$Expression_level)
# png(file="Plots/liver_3wa_E_gh_Slc11a1.png")
# Slc11a1.Egh_plot <- ggplot(Slc11a1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
#   geom_boxplot()+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25)+
#   labs(title = "Expression Level of Slc11a1 in Liver (E+B:gonad)", y = "Expression Level", x = "Genotype")+
#   theme_classic()+
#   ggeasy::easy_center_title()
# Slc11a1.Egh_plot
# dev.off()


Slc11a1.Egh <- g.h.E.reversed[,c("Slc11a1","genotype","treatment")]
colnames(Slc11a1.Egh)[1] <- "Expression_level"
Slc11a1.Egh$genotype <- as.factor(Slc11a1.Egh$genotype)
Slc11a1.Egh$Expression_level <- as.numeric(Slc11a1.Egh$Expression_level)

Slc11a1.Egh_plot <- ggplot(Slc11a1.Egh, aes(x= genotype, y= Expression_level, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Slc11a1", y = "Expression Level", x = "Genotype") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())+
  scale_fill_manual(values= colors_E,limits= c("E","B"))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5), dotsize = 0.4 ) 

## pairwise comparisons
Slc11a1.Egh.stat <- Slc11a1.Egh %>% 
  group_by(genotype) %>%
  t_test(Expression_level ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
Slc11a1.Egh.stat

Slc11a1.Egh.stat <- Slc11a1.Egh.stat %>%
  add_xy_position(x = "genotype")

Slc11a1.Egh_plot <- Slc11a1.Egh_plot+
  stat_pvalue_manual(
    Slc11a1.Egh.stat, label = "p.adj.signif", tip.lenghh = 0.01, hide.ns = TRUE, inherit.aes =FALSE)

ggsave("Plots_new/Liver_3wa_E_gh/Slc11a1.png", Slc11a1.Egh_plot, width = 5, height = 5, dpi = 300)


## arrange all plots

liver.3wa.E.gh.plots <- ggarrange(Slc11a1.Egh_plot,Sult3a1.Egh_plot, 
                                  labels = c("A", "B"),
                                  ncol = 2, 
                                  nrow = 1)
ggsave("Plots_new/Liver_3wa_E_gh/Combined.png", liver.3wa.E.gh.plots, width = 10, height = 5, dpi = 600)