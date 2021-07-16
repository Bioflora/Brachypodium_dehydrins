#Script_by_Pilar_Catalan_26/04/2021_MADR(13/07/2021)

# some of the listed packages below not necessary for PCA 
library("PerformanceAnalytics")
 library(ape)
 library(maps)
 library(phytools)
 library(diversitree)
 library(geiger)
 library(multcompView)
 library(caper)
 library(geiger)
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
 library(lmerTest)
 library(ggplot2)
 library(car)
 require(plotrix)
 library(caret)
 library(randomForest)
 library(effsize)
 library(ggstance)
 library(forcats)
library(viridis)
 library(Hmisc)


#import data 
climate<-read.csv(file="BdistachyonLongLatAltClimate.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
climate<-na.omit(climate)
clim2 = climate
clim2$ecotype<-as.factor(clim2$ecotype)

names(clim2)

#PCA analysis
pca_allBio <- dudi.pca(clim2[,5:23], scale = T, scan = FALSE, nf = 2)
screeplot(pca_allBio, main = "Screeplot - Eigenvalues")
vars<-get_pca_var(pca_allBio)
vars$contrib #variables contribution in the axis

#plotting ecotype contribution
fviz_pca_ind(pca_allBio, col.ind="contrib") +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=4)
#plotting variables contribution
fviz_pca_biplot(pca_allBio, title ="PCA climate variables", geom="point", pointsize = 2,  col.var="contrib")+
  scale_color_gradient2(low="deepskyblue1", mid="blue", high= "firebrick2", midpoint=2)+theme_grey()   #There were 50 or more warnings (use warnings() to see the first 50)


#summarize PCA1 climate by species
PCA1_sum <-aggregate(pca_allBio$li[,1], by=list(clim2$ecotype), FUN=mean, na.rm=TRUE)
PCA2_sum <-aggregate(pca_allBio$li[,2], by=list(clim2$ecotype), FUN=mean, na.rm=TRUE)


PCA1_sum
PCA2_sum

#extract PCA1 values
clim2$pca<-pca_allBio$li[,1]
pca_allBio$li[,1]


