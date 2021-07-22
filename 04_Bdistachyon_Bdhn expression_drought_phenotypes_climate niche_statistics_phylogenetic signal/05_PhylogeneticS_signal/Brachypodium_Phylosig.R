install.packages("ape")
install.packages("maps")
install.packages("phytools")
install.packages("geiger")
library(ape)
library(maps)
library(phytools)
library(geiger)

###############################################################################################################################################
######### Phylogenetic signal of Bdhn expression data (TPM), phenotypic changes and climate variation (PCA1) on BdistachyonBdhnTree #############
###############################################################################################################################################

###import Bdistachyon Bdhn tree and Bdhn_drought (D) data
BdistachyonBdhnTree.tree<-read.tree("BdistachyonBdhnTree.tre")
BdistachyonTPMD.data<-read.csv("BdistachyonTPM_D_BdhnTree.csv",row.names=1)
chk<-name.check(BdistachyonBdhnTree.tree,BdistachyonTPMD.data)
chk

#extract vectors (characters) from the data
Bdhn1aD <- as.matrix(BdistachyonTPMD.data)[, 1]
Bdhn2D <- as.matrix(BdistachyonTPMD.data)[, 2]
Bdhn3D <- as.matrix(BdistachyonTPMD.data)[, 3]
Bdhn7D <- as.matrix(BdistachyonTPMD.data)[, 4]

#check vectors
Bdhn1aD
Bdhn2D
Bdhn3D
Bdhn7D

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,BdistachyonTPMD.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPMD.data,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPMD.data,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))
lambda

###import Bdistachyon Bdhn tree and Bdhn_watered (W) data
BdistachyonBdhnTree.tree<-read.tree("BdistachyonBdhnTree.tre")
BdistachyonTPMW.data<-read.csv("BdistachyonTPM_W_BdhnTree.csv",row.names=1)
chk<-name.check(BdistachyonBdhnTree.tree,BdistachyonTPMW.data)
chk

#extract vectors (characters) from the data
Bdhn1aW <- as.matrix(BdistachyonTPMW.data)[, 1]
Bdhn2W <- as.matrix(BdistachyonTPMW.data)[, 2]
Bdhn3W <- as.matrix(BdistachyonTPMW.data)[, 3]
Bdhn7W <- as.matrix(BdistachyonTPMW.data)[, 4]

#check vectors
Bdhn1aW
Bdhn2W
Bdhn3W
Bdhn7W

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,BdistachyonTPMW.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPMW.data,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPMW.data,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))
lambda

###import Bdistachyon Bdhn tree and phenotypic traits_watered (W) & drought (D) data and climate PCA1 data

BdistachyonBdhnTree.tree<-read.tree("BdistachyonBdhnTree.tre")
BdistachyonPhenoClimBdhnTree.data<-read.csv("BrachyPhenoW&D_ClimPCA1_BdhnTree.csv",row.names=1)
chk<-name.check(BdistachyonBdhnTree.tree,BdistachyonPhenoClimBdhnTree.data)
chk

#extract vectors (characters) from the data
leaf_rwcW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 1]
leaf_wcW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 2]
lmaW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 3]
proW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 4]
abvrgdW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 5]
blwgrdW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 6]
ttlmassW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 7]
rmrW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 8]
delta13c2W <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 9]
leafcW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 10]
leafnW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 11]
cnW <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 12]
leaf_rwcD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 13]
leaf_wcD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 14]
lmaD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 15]
proD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 16]
abvrgdD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 17]
blwgrdD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 18]
ttlmassD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 19]
rmrD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 20]
delta13cD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 21]
leafcD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 22]
leafnD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 23]
cnD <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 24]
climPCA1 <- as.matrix(BdistachyonPhenoClimBdhnTree.data)[, 25]

#check vectors
leaf_rwcW
leaf_wcW
lmaW
proW
abvrgdW
blwgrdW
ttlmassW
rmrW
delta13c2W
leafcW
leafnW
cnW
leaf_rwcD
leaf_wcD
lmaD
proD
abvrgdD
blwgrdD
ttlmassD
rmrD
delta13cD
leafcD
leafnD
cnD
climPCA1

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonBdhntree.u<-chronopl(BdistachyonBdhnTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonBdhntree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonBdhntree.u,BdistachyonPhenoClimBdhnTree.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonPhenoClimBdhnTree.data,2,phylosig,tree=BdistachyonBdhntree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonPhenoClimBdhnTree.data,2,phylosig,tree=BdistachyonBdhntree.u,method="lambda",test=TRUE)))
lambda

###############################################################################################################################################
######### Phylogenetic signal of Bdhn expression data (TPM), phenotypic changes and climate variation (PCA1) on BdistachyonSpTree #############
###############################################################################################################################################

###import Bdistachyon species tree and Bdhn_drought (D) data
BdistachyonSpTree.tree<-read.tree("BdistachyonSpTree.tre")
BdistachyonTPMD.data<-read.csv("BdistachyonTPM_D_SpTree.csv",row.names=1)
chk<-name.check(BdistachyonSpTree.tree,BdistachyonTPMD.data)
chk

#extract vectors (characters) from the data
Bdhn1aD <- as.matrix(BdistachyonTPMD.data)[, 1]
Bdhn2D <- as.matrix(BdistachyonTPMD.data)[, 2]
Bdhn3D <- as.matrix(BdistachyonTPMD.data)[, 3]
Bdhn7D <- as.matrix(BdistachyonTPMD.data)[, 4]

#check vectors
Bdhn1aD
Bdhn2D
Bdhn3D
Bdhn7D

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonSptree.u<-chronopl(BdistachyonSpTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonSptree.u,BdistachyonTPMD.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPMD.data,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPMD.data,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))
lambda

###import Bdistachyon species tree and Bdhn_watered (W) data
> BdistachyonSpTree.tree<-read.tree("BdistachyonSpTree.tre")
> BdistachyonTPMW.data<-read.csv("BdistachyonTPM_W_SpTree.csv",row.names=1)
> chk<-name.check(BdistachyonSpTree.tree,BdistachyonTPMW.data)
> chk

#extract vectors (characters) from the data
Bdhn1aW <- as.matrix(BdistachyonTPMW.data)[, 1]
Bdhn2W <- as.matrix(BdistachyonTPMW.data)[, 2]
Bdhn3W <- as.matrix(BdistachyonTPMW.data)[, 3]
Bdhn7W <- as.matrix(BdistachyonTPMW.data)[, 4]

#check vectors
Bdhn1aW
Bdhn2W
Bdhn3W
Bdhn7W

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonSptree.u<-chronopl(BdistachyonSpTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonSptree.u,BdistachyonTPMW.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonTPMW.data,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonTPMW.data,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))
lambda

###import Bdistachyon species tree and phenotypic traits_watered (W) & drought (D) data and climate PCA1 data

BdistachyonSpTree.tree<-read.tree("BdistachyonSpTree.tre")
BdistachyonPhenoClimSpTree.data<-read.csv("BrachyPhenoW&D_ClimPCA1_SpTree.csv",row.names=1)
chk<-name.check(BdistachyonSpTree.tree,BdistachyonPhenoClimSpTree.data)
chk

#extract vectors (characters) from the data
leaf_rwcW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 1]
leaf_wcW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 2]
lmaW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 3]
proW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 4]
abvrgdW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 5]
blwgrdW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 6]
ttlmassW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 7]
rmrW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 8]
delta13c2W <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 9]
leafcW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 10]
leafnW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 11]
cnW <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 12]
leaf_rwcD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 13]
leaf_wcD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 14]
lmaD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 15]
proD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 16]
abvrgdD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 17]
blwgrdD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 18]
ttlmassD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 19]
rmrD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 20]
delta13cD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 21]
leafcD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 22]
leafnD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 23]
cnD <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 24]
climPCA1 <- as.matrix(BdistachyonPhenoClimSpTree.data)[, 25]

#check vectors
leaf_rwcW
leaf_wcW
lmaW
proW
abvrgdW
blwgrdW
ttlmassW
rmrW
delta13c2W
leafcW
leafnW
cnW
leaf_rwcD
leaf_wcD
lmaD
proD
abvrgdD
blwgrdD
ttlmassD
rmrD
delta13cD
leafcD
leafnD
cnD
climPCA1

#ultrametric tree (branches could not have zero lenght) (age-max = 1 Ma, based on Sancho et al. 2018)
BdistachyonSptree.u<-chronopl(BdistachyonSpTree.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
plot(BdistachyonSptree.u)

#phylo.heatmap
phylo.heatmap(BdistachyonSptree.u,BdistachyonPhenoClimSpTree.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#phylogenetic signal, phylosig_K
K<-apply(BdistachyonPhenoClimSpTree.data,2,phylosig,tree=BdistachyonSptree.u, method="K", test=TRUE, nsim=1000)
K

#phylogenetic signal, phylosig_lambda
lambda<-t(simplify2array(apply(BdistachyonPhenoClimSpTree.data,2,phylosig,tree=BdistachyonSptree.u,method="lambda",test=TRUE)))
lambda

