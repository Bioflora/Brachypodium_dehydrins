install.packages("ape")
install.packages("maps")
install.packages("phytools")
install.packages("geiger")
library(ape)
library(maps)
library(phytools)
library(geiger)


###############################################################################################################################################
######### Phylogenetic signal of Bdhn expression data (TPM), phenotypic changes and climate variation (PCA1) on BdistachyonBdhnTree ###########
###############################################################################################################################################
##Fig 6a.
###import Bdistachyon Bdhn tree and Bdhn_drought (D) data
BdistachyonBdhnTree.tree<-read.tree("BdistachyonBdhnTree.tre")
BdistachyonTPMD.data<-read.csv("BdistachyonTPM_D_BdhnTree.csv",row.names=1)
BdistachyonTPMW.data<-read.csv("BdistachyonTPM_W_BdhnTree.csv",row.names=1)
BdistachyonTPM = cbind( BdistachyonTPMW.data, BdistachyonTPMD.data)
chk<-name.check(BdistachyonBdhnTree.tree,BdistachyonTPM)


#extract vectors (characters) from the data
Bdhn1aW <- as.matrix(BdistachyonTPM)[, 1]
Bdhn2W <- as.matrix(BdistachyonTPM)[, 2]
Bdhn3W <- as.matrix(BdistachyonTPM)[, 3]
Bdhn7W <- as.matrix(BdistachyonTPM)[, 4]
Bdhn1aD <- as.matrix(BdistachyonTPM)[, 5]
Bdhn2D <- as.matrix(BdistachyonTPM)[, 6]
Bdhn3D <- as.matrix(BdistachyonTPM)[, 7]
Bdhn7D <- as.matrix(BdistachyonTPM)[, 8]

#check vectors
Bdhn1aD
Bdhn1aW

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonBdhntree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,BdistachyonTPM,standardize=TRUE,fsize=c(0.9,0.9,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPM,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPM,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))

###import Bdistachyon Bdhn tree and phenotypic traits_watered (W) & drought (D) data and climate PCA1 data
## Fig. 6b

BdistachyonBdhnTree.tree<-read.tree("BdistachyonBdhnTree.tre")
BdistachyonPhenoClimBdhnTree.data<-read.csv("BrachyPhenoW&D_ClimPCA1_BdhnTree.csv",row.names=1)
name.check(BdistachyonBdhnTree.tree,BdistachyonPhenoClimBdhnTree.data)
names(BdistachyonPhenoClimBdhnTree.data)[5] <- "abvgrdW"
names(BdistachyonPhenoClimBdhnTree.data)[17] <- "abvgrdD"

##Split dataset into to datasets
watered<- subset(BdistachyonPhenoClimBdhnTree.data, select=c(1:12))
drought<- subset(BdistachyonPhenoClimBdhnTree.data, select=c(13:24))

#extract vectors (characters) from the data
leaf_rwcW <- as.matrix(watered)[, 1]
leaf_wcW <- as.matrix(watered)[, 2]
lmaW <- as.matrix(watered)[, 3]
proW <- as.matrix(watered)[, 4]
abvrgdW <- as.matrix(watered)[, 5]
blwgrdW <- as.matrix(watered)[, 6]
ttlmassW <- as.matrix(watered)[, 7]
rmrW <- as.matrix(watered)[, 8]
delta13c2W <- as.matrix(watered)[, 9]
leafcW <- as.matrix(watered)[, 10]
leafnW <- as.matrix(watered)[, 11]
cnW <- as.matrix(watered)[, 12]
leaf_rwcD <- as.matrix(drought)[, 1]
leaf_wcD <- as.matrix(drought)[, 2]
lmaD <- as.matrix(drought)[, 3]
proD <- as.matrix(drought)[, 4]
abvrgdD <- as.matrix(drought)[, 5]
blwgrdD <- as.matrix(drought)[, 6]
ttlmassD <- as.matrix(drought)[, 7]
rmrD <- as.matrix(drought)[, 8]
delta13cD <- as.matrix(drought)[, 9]
leafcD <- as.matrix(drought)[, 10]
leafnD <- as.matrix(drought)[, 11]
cnD <- as.matrix(drought)[, 12]


#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonBdhntree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,watered,standardize=TRUE,fsize=c(0.9,0.9,1))
phylo.heatmap(BdistachyonBdhntree.u,drought,standardize=TRUE,fsize=c(0.9,0.9,1))

#phylogenetic signal, phylosig_K
Kw<-apply(watered,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)
Kd<-apply(drought,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)

#phylogenetic signal, phylosig_lambda
lambdaw<-t(simplify2array(apply(watered,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))
lambdad<-t(simplify2array(apply(drought,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))


####PCA1######
#fig.6c
PCA<-read.csv("BdistachyonClimPCA1mod_BdhnTree.csv",row.names=1)
names(PCA)[1] <- "PCA1"
name.check(BdistachyonBdhnTree.tree,PCA)
PCA1 <- as.matrix(PCA)[, 1]
aux <- as.matrix(PCA)[, 2]

#ultrametric tree 
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonBdhntree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,PCA,standardize=TRUE,fsize=c(0.9,0.9,1))

###############################################################################################################################################
######### Phylogenetic signal of Bdhn expression data (TPM), phenotypic changes and climate variation (PCA1) on BdistachyonSpTree #############
###############################################################################################################################################

###import Bdistachyon species tree and Bdhn_drought (D) data
BdistachyonSpTree.tree<-read.tree("BdistachyonSpTree.tre")
BdistachyonTPMD.data<-read.csv("BdistachyonTPM_D_SpTree.csv",row.names=1)
BdistachyonTPMW.data<-read.csv("BdistachyonTPM_W_SpTree.csv",row.names=1)
BdistachyonTPM_Bdhn_SP<- cbind(BdistachyonTPMW.data, BdistachyonTPMD.data )
name.check(BdistachyonSpTree.tree,BdistachyonTPM_Bdhn_SP)


#extract vectors (characters) from the data
Bdhn1aW <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 1]
Bdhn2W <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 2]
Bdhn3W <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 3]
Bdhn7W <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 4]
Bdhn1aD <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 5]
Bdhn2D <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 6]
Bdhn3D <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 7]
Bdhn7D <- as.matrix(BdistachyonTPM_Bdhn_SP)[, 8]
#check vectors

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonSptree.u<-chronopl(BdistachyonSpTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonSptree.u,BdistachyonTPM_Bdhn_SP,standardize=TRUE,fsize=c(0.9,0.9,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPM_Bdhn_SP,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)


#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPM_Bdhn_SP,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))


################################################################################################################
###import Bdistachyon species tree and phenotypic traits_watered (W) & drought (D) data and climate PCA1 data
##Supplementary figures

BdistachyonSpTree.tree<-read.tree("BdistachyonSpTree.tre")
BdistachyonPhenoClimSpTree.data<-read.csv("BrachyPhenoW&D_ClimPCA1_SpTree.csv",row.names=1)
PCA1<- read.csv("ClimPCA1_SpTree2.csv", row.names=1)
name.check(BdistachyonSpTree.tree,PCA1)
watered_sp<- subset(BdistachyonPhenoClimSpTree.data, select=c(1:12))
drought_sp<- subset(BdistachyonPhenoClimSpTree.data, select=c(13:24))

#extract vectors (characters) from the data
leaf_rwcW <- as.matrix(watered_sp)[, 1]
leaf_wcW <- as.matrix(watered_sp)[, 2]
lmaW <- as.matrix(watered_sp)[, 3]
proW <- as.matrix(watered_sp)[, 4]
abvrgdW <- as.matrix(watered_sp)[, 5]
blwgrdW <- as.matrix(watered_sp)[, 6]
ttlmassW <- as.matrix(watered_sp)[, 7]
rmrW <- as.matrix(watered_sp)[, 8]
delta13c2W <- as.matrix(watered_sp)[, 9]
leafcW <- as.matrix(watered_sp)[, 10]
leafnW <- as.matrix(watered_sp)[, 11]
cnW <- as.matrix(watered_sp)[, 12]
#Drought
leaf_rwcD <- as.matrix(drought_sp)[, 1]
leaf_wcD <- as.matrix(drought_sp)[, 2]
lmaD <- as.matrix(drought_sp)[, 3]
proD <- as.matrix(drought_sp)[, 4]
abvrgdD <- as.matrix(drought_sp)[, 5]
blwgrdD <- as.matrix(drought_sp)[, 6]
ttlmassD <- as.matrix(drought_sp)[, 7]
rmrD <- as.matrix(drought_sp)[, 8]
delta13cD <- as.matrix(drought_sp)[, 9]
leafcD <- as.matrix(drought_sp)[, 10]
leafnD <- as.matrix(drought_sp)[, 11]
cnD <- as.matrix(drought_sp)[, 12]
#PCA1
pca1<- as.matrix(PCA1) [, 1]

#check vectors

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonSptree.u<-chronopl(BdistachyonSpTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap_W
phylo.heatmap(BdistachyonSptree.u,watered_sp,standardize=TRUE,fsize=c(0.9,0.9,1))
phylo.heatmap(BdistachyonSptree.u,drought_sp,standardize=TRUE,fsize=c(0.9,0.9,1))
phylo.heatmap(BdistachyonSptree.u,PCA1,standardize=TRUE,fsize=c(0.9,0.9,1))

#phylogenetic signal, phylosig_K
Kw<-apply(watered_sp,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)
Kd<-apply(drought_sp,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)
Kpca<-apply(PCA1,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000) 

#phylogenetic signal, phylosig_lambda
lambdaw<-t(simplify2array(apply(watered_sp,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))
lambdad<-t(simplify2array(apply(drought_sp,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))
lambdapca<-t(simplify2array(apply(PCA1,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE))) 
