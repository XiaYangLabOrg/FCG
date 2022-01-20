library("genefilter")
ing <- as.matrix(read.table("expr.lumi.gene(8751).log2.0.01.dr0.5.txt", header = T))
ing.var <- varFilter(ing, var.cutoff = 0.6)
ing.cor <- ing.var[,c("X5719882029_E", "X5613711027_B", "X5719882029_F", "X5613711027_C", "X5613711015_H", "X5773660159_G")]
cor <- cor(ing.cor)
plot(density(cor[cor!=1]), main = "10539 in XXSry - 40%")
abline(v = cor["X5613711027_B", "X5719882029_E"])
ing.cor.xy <- ing.var[,c("X5719882029_E", "X5613711027_B", "X5613711015_G", "X5773660104_E", "X5773660159_F")]
cor.xy <- cor(ing.cor.xy)
plot(density(cor.xy[cor.xy!=1]), main = "10539 in XY-Sry - 40%")
abline(v = cor["X5613711027_B", "X5719882029_E"])

XXsry.cor <- ing.var[,c("X5773660094_B", "X5773660104_B", "X5613711014_E", "X5773660134_C", "X5613711027_G")]
XXsry <- cor(XXsry.cor)
plot(density(XXsry[XXsry!=1]), main = "75533 - 40%")
abline(v = XXsry["X5613711014_E", "X5773660134_C"])

XX.cor <- ing.var[,c("X5773660159_A", "X5719882029_H", "X5773660104_H", "X5613711015_B")]
XX <- cor(XX.cor)
plot(density(XX[XX!=1]), main = "69323 - 40%")
abline(v = XX["X5719882029_H", "X5773660104_H"])

liver <- as.matrix(read.table("../FCG/liver/expr.lumi.gene(6989).log2.0.01.dr0.5.txt", header = T))
liver.var <- varFilter(liver, var.cutoff = 0.6)
liver.cor <- liver.var[,c("X5613703028_F", "X5719882028_C", "X5719882013_A", "X5719882030_H", "X5719882006_D")]
cor.liver <- cor(liver.cor)
