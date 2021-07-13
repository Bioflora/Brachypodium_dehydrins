#Script_by_Pilar_Catalan_26/04/2021_MADR(5/05/2021)


# some of the listed packages below not necessary for PCA 

install.packages("ape")
install.packages("maps")
install.packages("caper")
install.packages("phytools")
install.packages("diversitree")
install.packages("geiger")
install.packages("multcompView")
install.packages("lme4")
install.packages("nlme")
install.packages("vegan")
install.packages("factoextra")
install.packages("ade4")
install.packages("ggfortify")
install.packages("corrplot")
install.packages("lmerTest")
install.packages("lsmeans")
install.packages("agricolae")
install.packages("multcomp")
install.packages("pbkrtest")
install.packages("ggplot2")
install.packages("car")
install.packages("plotrix")
install.packages("caret")
install.packages("randomForest")
install.packages("effsize")
install.packages("ggstance")
install.packages("forcats")
install.packages("Hmisc")
install.packages("MCMCglmm")
install.packages("extRemes")
install.packages("MuMIn")
install.packages("adephylo")
install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
 library(ape)
 library(maps)
 library(phytools)
 library(diversitree)
 library(geiger)
 library(multcompView)
 library(caper)
 library(caper)
 library(phytools)
 library(diversitree)
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
library(viridis)
 library(Hmisc)


#import data 
climate<-read.csv(file="BdisLongLatAlt.csv", h=T, sep=";", na.strings = "NA", stringsAsFactors=FALSE)
climate<-na.omit(climate)
clim2 = climate
clim2$ecotype<-as.factor(clim2$ecotype)

names(clim2)
clim2
chart.Correlation(clim2 [,6:24], histogram=TRUE, pch=19)

#PCA analysis
pca_allBio <- dudi.pca(clim2[,6:24], scale = T, scan = FALSE, nf = 2)
screeplot(pca_allBio, main = "Screeplot - Eigenvalues")
vars<-get_pca_var(pca_allBio)

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

#summarize PCA1 climate by species
PCA1_sum <-aggregate(pca_allBio$li[,1], by=list(clim2$ecotype), FUN=mean, na.rm=TRUE)
PCA1_sum


