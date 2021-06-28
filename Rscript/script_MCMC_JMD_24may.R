require(caper)
require(phytools)
require(diversitree)
library(phylotools)
library(ape)
library(geiger)
library(multcompView)
library(lme4)
library(nlme)
library(vegan)
library(factoextra)
library(ade4)
library(ggfortify)
library(corrplot)
library(lmerTest)
library(lsmeans)
library(agricolae)
library(multcomp)
library(pbkrtest)
library(vegan)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(nlme)
library(car)
require(plotrix)
library(caret)
library(randomForest)
library(effsize)
library(ggstance)
library(forcats)
library(Hmisc)
library(MCMCglmm)
library(extRemes)
library( MuMIn )
library(ggplot2)
library(adephylo)
library(effsize)
library(tidyverse)
library(maps)
library(dplyr)
library(readxl)
library(ggsignif)
library(ggpubr)

# import data

tab1<-read.csv(file="BrachyBdhnPhenoW&DallmodABR8noBdTR3cnoBdTR13a.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
tab2<- na.omit(tab_all)
delta13c= as.data.frame(tab2 [19, drop=FALSE]) #extrae la columna
delta13C<- as.numeric(tab2$delta13c)
delta13c_abs<-abs(delta13C) #hace valores absolutos
aux1 <- subset(tab2, select=-c(19)) #quita la columna
tab3<- cbind(aux1, delta13c_abs) #une la columna mod a la anterior
tab4<- slice(tab3, - (38:44)) #elimina ABR8
tab4$ecotype<-as.factor(tab4$ecotype)
tab4$treatment<-as.factor(tab4$treatment)
tab4$climate<-as.factor(tab4$climate)
names(tab4)

#Transforming into Ultrametric tree (branches could not have zero lenght. )
BdhnTreeD1<-read.tree(file="BrachyBdhnTreeDMD2.tree") #Branches with 0 lenght were slightly modify with 0.000000001, manually)
plot(BdhnTreeD1)
pr.BdhntreeD.u<-chronopl(BdhnTreeD1, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(pr.BdhntreeD.u)
Bdhn.u <-pr.BdhntreeD.u

## compute Hedges' g Bdhn and other variables###
tab4$animal<-tab4$ecotype
tab4
tab4$animal<-as.factor(tab4$animal)
tab4$animal
levels1<- as.factor(levels(tab4$animal))
levels1

#Matriz en blanco (ajusta a 29 despues de quitar ABR8)
mat_rb<-matrix(rep(NA,29*3),ncol=3,nrow=29); colnames(mat_rb)= c("mean", "inf", "sup"); rownames(mat_rb)= levels(tab4$animal)
mat_rb1<- mat_rb

##cohen.d analysis #Matriz para cada analisis guardada con su nombre para usarla posteriormente

for (i in 1:29){
  es<-cohen.d(subset(tab4, animal == levels(tab4$animal)[i])$delta13c_abs,
              subset(tab4, animal == levels(tab4$animal)[i])$treatment,hedges.correction=TRUE, pooled = T)      
  mat_rb1[i,1]<-es$estimate
  mat_rb1[i,2]<-es$conf.int[1]
  mat_rb1[i,3]<-es$conf.int[2]
}
mat1<-mat_rb1 
mean(mat1[,1])
mean(mat2[,1])
mean(mat3[,1])
mean(mat7[,1])



#DOTCHART
#######Para Bdhn 
df1<-as.data.frame(mat1)
df1$ecotype<-rownames(mat1)
df2<-df1[match(Bdhn.u$tip.label, df1$ecotype), ] #arbol ultrametrico  
plot(Bdhn.u)



title.1 <- "Effect of delta13c with 95% c.i."
title.2 <- "Effect of leaf_rwc with 95% c.i."
df2= -1*df2#no entiendo que hace esto porque lo que consigue es eliminar el ecotipo que antes tenia

dotchart(df2$mean, col= c("#1E90FF"), 
         pch = 19, cex = 1, xlim = c(-25,2), labels=row.names(df2), main = title.1, cex.main= 1)
abline(v=0, col="black", lty = 2)
abline(v=mean(df2[,1]), col ="red", lty=2)

for (i in 1:nrow(df2)){
  lines(x=c(df2$inf[i],df2$sup[i]), y=c(i,i))
  col = "grey"
}
grid()




####  MCMC   #######
####################
names(tab4)
tab4$rb<-log(tab4$lma)
tab5<- as.data.frame(tab4)
Bdhn.u$tip.label
Bdhn.u<-makeLabel(Bdhn.u)
#plot(Bdhn.u)

#0Bdhn7$animal<-as.factor(Bdhn7$animal)
#Bdhn7$animal
#levels<- as.factor(levels(Bdhn7$animal))
#levels
##
summary(test1<-MCMCglmm(rb~treatment , random=~animal, pedigree = Bdhn.u,  data =tab4,   family = "gaussian", singular.ok=TRUE))

boxplot(tab4$rb~tab4$treatment, data = tab4, col = c("firebrick1", "cyan3"), xlab= "treatment", ylab = "gm-2", main="lma", cex.lab= 1.5, cex.axis=1.3); text(1.5,8000, "*", cex = 3)

###########################
# Load data
leaf_rwc<- subset.data.frame(tab4, select=c(4,11))
leaf_wc<- subset.data.frame(tab4, select=c(4,12))
lma<- subset.data.frame(tab4, select=c(4,13))
pro<- subset.data.frame(tab4, select=c(4,14))
abvgrd<- subset.data.frame(tab4, select=c(4,15))
blwgrd<- subset.data.frame(tab4, select=c(4,16))
ttlmass<- subset.data.frame(tab4, select=c(4,17))
rmr<- subset.data.frame(tab4, select=c(4,18))
leafc<- subset.data.frame(tab4, select=c(4,19))
leafn<- subset.data.frame(tab4, select=c(4,20))
cn<- subset.data.frame(tab4, select=c(4,21))
WUE<- subset.data.frame(tab4, select=c(4,22))


# Plot in png

png("All_phenotraits_DW.png", width=1000, height=800)

leaf_rwc_DW <- ggplot(data = leaf_rwc, aes(x=treatment, y=leaf_rwc)) + geom_boxplot(aes(fill=treatment)) + ggtitle("Leaf_rwc") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(120), map_signif_level=TRUE)

leaf_wc_DW <- ggplot(data = leaf_wc, aes(x=treatment, y=leaf_wc)) + geom_boxplot(aes(fill=treatment)) + ggtitle("leaf_wc") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(400), map_signif_level=TRUE)


lma_DW <- ggplot(data = lma, aes(x=treatment, y=lma)) + geom_boxplot(aes(fill=treatment)) + ggtitle("lma") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(40), map_signif_level=TRUE)

pro_DW <- ggplot(data = pro, aes(x=treatment, y=pro)) + geom_boxplot(aes(fill=treatment)) + ggtitle("pro") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(100), map_signif_level=TRUE)

Abvgrd_DW <- ggplot(data = abvgrd, aes(x=treatment, y=abvrgd)) + geom_boxplot(aes(fill=treatment)) + ggtitle("abvgrd") + theme(plot.title = element_text(hjust = 0.5)) +
geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(250), map_signif_level=TRUE)

blwgrd_DW <- ggplot(data = blwgrd, aes(x=treatment, y=blwgrd)) + geom_boxplot(aes(fill=treatment)) + ggtitle("blwgrd") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(100), map_signif_level=TRUE)

ttlmass_DW <- ggplot(data = ttlmass, aes(x=treatment, y=ttlmass)) + geom_boxplot(aes(fill=treatment)) + ggtitle("ttlmass") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(350), map_signif_level=TRUE)

rmr_DW <- ggplot(data = rmr, aes(x=treatment, y=rmr)) + geom_boxplot(aes(fill=treatment)) + ggtitle("rmr") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(50), map_signif_level=TRUE)

WUE_DW <- ggplot(data = WUE, aes(x=treatment, y=delta13c_abs)) + geom_boxplot(aes(fill=treatment)) + ggtitle("WUE") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(35), map_signif_level=TRUE)

leafc_DW <- ggplot(data = leafc, aes(x=treatment, y=leafc)) + geom_boxplot(aes(fill=treatment)) + ggtitle("leafc") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(420), map_signif_level=TRUE)

leafn_DW <- ggplot(data =leafn, aes(x=treatment, y=leafn)) + geom_boxplot(aes(fill=treatment)) + ggtitle("leafn") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(50), map_signif_level=TRUE)

cn_DW <- ggplot(data = cn, aes(x=treatment, y=cn)) + geom_boxplot(aes(fill=treatment)) + ggtitle("C:N") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Drought","Watered")), y_position = c(20), map_signif_level=TRUE)



ggarrange(leaf_rwc_DW, leaf_wc_DW, lma_DW, pro_DW, Abvgrd_DW, blwgrd_DW, ttlmass_DW, rmr_DW, WUE_DW, leafc_DW, leafn_DW, cn_DW,  ncol = 4, nrow = 3, common.legend = TRUE, legend="bottom")


dev.off()


