
#Script_by_Pilar_Catalan_26/04/2021

> install.packages("ape")
> install.packages("maps")
> install.packages("phytools")
> library(ape)
> library(maps)
> library(phytools)
> library(geiger)

#import tree and data
> BrachySpTreeD.tree<-read.tree("BrachySpTreeD.tre")
> BrachyTPMD.data<-read.csv("BrachyTPM_D_SpTree.csv",row.names=1)
> chk<-name.check(BrachySpTreeD.tree,BrachyTPMD.data)
> chk
[1] "OK"

#extract a vector (character) from the data
Bdhn7 <- as.matrix(BrachyTPMD.data)[, 4]
Bdhn7
       ABR2        ABR3        ABR4        ABR5        ABR6 
       19.4        29.5        15.1        44.2        41.6 
      Adi10       Adi12        Adi2       Bd1-1      Bd18-1 
      332.9       139.8        63.8        25.6       172.9 
Bd21control      Bd21-3       Bd2-3      Bd30-1       Bd3-1 
       24.3       296.5        29.3       222.7       329.4 
    BdTR10c     BdTR11g     BdTR11i     BdTR12c      BdTR1i 
      194.1        95.7       124.3         5.8       251.2 
     BdTR2b      BdTR2g      BdTR5i      BdTR9k        Bis1 
      136.4       517.7        70.7        78.4       149.6 
       Kah1        Kah5        Koz1        Koz3        RON2 
       65.9        86.4        57.4       812.2       140.4 

#Ultrametric tree (branches could not have zero lenght)
> pr.tree.u<-chronopl(BrachyBdhnTreeD.tree, lambda=0.5, age.max = 1,node = "root", tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
> plot(pr.tree.u)

#Phylo.heatmap
> phylo.heatmap(BrachySpTreeD.tree,BrachyTPMD.data,standardize=TRUE,fsize=c(0.5,0.7,1))

#Phylosig_K_&_lambda_datamatrix
> K<-apply(BrachyTPMD.data,2,phylosig,tree=BrachySpTreeD.tree)
> K
    Bdhn1aD      Bdhn2D      Bdhn3D      Bdhn7D 
0.003377399 0.003063527 0.002240404 0.001131493 

>K<-apply(BrachyPhenoD.data,2,phylosig,tree=pr.tree.u, method="K", test=TRUE, nsim=1000)
> K
$leaf_rwcD
Phylogenetic signal K : 0.00595303 
P-value (based on 1000 randomizations) : 0.186 

$leaf_wcD
Phylogenetic signal K : 0.0191719 
P-value (based on 1000 randomizations) : 0.035 

$lmaD
Phylogenetic signal K : 0.00655413 
P-value (based on 1000 randomizations) : 0.192
 
..............

> lambda<-t(simplify2array(apply(BrachyTPMD.data,2,phylosig,tree=BrachySpTreeD.tree,method="lambda")))
> lambda
        lambda       logL      lik
Bdhn1aD 0.5671077    -194.0138 ?  
Bdhn2D  4.526224e-05 -211.6308 ?  
Bdhn3D  6.610696e-05 -238.9317 ?  
Bdhn7D  6.610696e-05 -197.27   ?

> lambda<-t(simplify2array(apply(BrachyPhenoD.data,2,phylosig,tree=pr.tree.u,method="lambda",test=TRUE)))
> lambda
          lambda       logL      logL0     P           lik
leaf_rwcD 0.8240433    -70.83287 -71.94764 0.1353942   ?  
leaf_wcD  0.5318181    -129.0792 -129.2025 0.6194317   ?  
lmaD      0.4014225    -60.3207  -60.04825 1           ?  
proD      6.611279e-05 -118.2053 -118.2048 1           ?  
abvrgdD   0.9862024    -132.3802 -137.5324 0.001327115 ?  
blwgrdD   0.6778591    -112.6454 -113.232  0.2787568   ?  
ttlmassD  0.9586421    -144.0437 -147.1176 0.01315614  ?  
rmrD      0.8587116    -71.60054 -74.00027 0.02846858  ?  
delta13cD 0.5455282    -9.438273 -9.102205 1           ?  
leafcD    6.611279e-05 -86.81046 -86.80988 1           ?  
leafnD    0.8379002    -52.22653 -57.04403 0.001909058 ?  
cnD       0.8350859    -24.49639 -28.06983 0.007509447 ?  


#Phylosig_K_one_character
> phylosig(pr.BrachySptree.u,Bdhn7,method="K", test=TRUE, nsim=1000)
Phylogenetic signal K : 0.0476759 
P-value (based on 1000 randomizations) : 0.327   #no phylogenetic signal

> phylosig(pr.BrachySptree.u,BdhnF,method="K", test=TRUE, nsim=1000)
Phylogenetic signal K : 6.22598 
P-value (based on 1000 randomizations) : 0.001   #phylogenetic signal (p=0.001)

#Phylosig_lambda_one_character
> phylosig(pr.BrachySptree.u,Bdhn7,method="lambda", test=TRUE)
Phylogenetic signal lambda : 6.61222e-05 
logL(lambda) : -196.53 
LR(lambda=0) : -0.00105877 
P-value (based on LR test) : 1     #no phylogenetic signal 

> phylosig(pr.BrachySptree.u,BdhnF,method="lambda", test=TRUE)
Phylogenetic signal lambda : 1.00016 
logL(lambda) : -109.32 
LR(lambda=0) : 65.2108 
P-value (based on LR test) : 6.72983e-16 #phylogenetic signal

#Ancestral States 
> plotTree(BrachySpTreeD.tree)
> BrachyBdhn1aDSpTree.data<-read.csv("BrachyBdhn1a_D_SpTree.csv",row.names=1)
> BrachyBdhn1aDSpTree<-read.csv("BrachyBdhn1a_D_SpTree.csv",row.names=1)
> BrachyBdhn1aD<-as.matrix(BrachyBdhn1aDSpTree)[,1]
> BrachyBdhn1aD
       ABR2        ABR3        ABR4        ABR5        ABR6       Adi10       Adi12 
      214.2       263.8       247.2       243.9       290.2       659.4       488.7 
       Adi2       Bd1-1      Bd18-1 Bd21control      Bd21-3       Bd2-3      Bd30-1 
      460.0       403.1       737.0       326.9       600.4       502.2       488.7 
      Bd3-1     BdTR10c     BdTR11g     BdTR11i     BdTR12c      BdTR1i      BdTR2b 
      395.2       410.4       497.2       531.9       297.8       582.9       538.2 
     BdTR2g      BdTR5i      BdTR9k        Bis1        Kah1        Kah5        Koz1 
      750.6       359.6       436.8       473.9       385.4       432.8       419.6 
       Koz3        RON2 
      898.5       656.2 
> fit<-fastAnc(BrachySpTreeD.tree,BrachyBdhn1aD,vars=TRUE,CI=TRUE)
> fit
Ancestral character estimates using fastAnc:
      31       32       33       34       35       36       37       38       39       40 
391.5653 390.3518 458.6395 449.0911 466.7612 467.5705 478.5956 429.1232 445.0457 461.3959 
      41       42       43       44       45       46       47       48       49       50 
459.2563 442.7445 466.5386 485.5514 480.6111 513.2337 517.9404 563.1807 603.9712 581.0390 
      51       52       53       54       55       56       57       58       59 
620.2997 635.5329 592.5777 417.2185 373.9223 357.2034 357.6604 373.4462 334.2881 

Variances on ancestral states:
          31           32           33           34           35           36           37 
1736832.0960  792469.5440  743728.3626  700212.6633  707003.8178  478569.9684  535245.5795 
          38           39           40           41           42           43           44 
 487975.1427  425326.8464  271387.0612  431139.3468  488752.5006  452885.1343  436839.1421 
          45           46           47           48           49           50           51 
 546134.0212     445.6194  468914.1269  512747.2086  630204.1831  611556.8636     780.2528 
          52           53           54           55           56           57           58 
    364.5537  753194.4437  903088.9236  794717.7969  741478.3593  739208.0804  722502.2161 
          59 
 709950.5145 

Lower & upper 95% CIs:
        lower     upper
31 -2191.4976 2974.6282
32 -1354.4550 2135.1587
33 -1231.6584 2148.9373
34 -1191.0116 2089.1938
35 -1181.2758 2114.7981
36  -888.3331 1823.4740
37  -955.3498 1912.5410
38  -940.0391 1798.2854
39  -833.2092 1723.3006
40  -559.6626 1482.4545
41  -827.7033 1746.2159
42  -927.5079 1812.9969
43  -852.4776 1785.5547
44  -809.8872 1780.9901
45  -967.8462 1929.0684
46   471.8587  554.6087
47  -824.2148 1860.0956
48  -840.3041 1966.6656
49  -951.9826 2159.9251
50  -951.7222 2113.8002
51   565.5510  675.0484
52   598.1100  672.9557
53 -1108.4431 2293.5985
54 -1445.3889 2279.8259
55 -1373.3578 2121.2024
56 -1330.5357 2044.9425
57 -1327.4930 2042.8137
58 -1292.5563 2039.4488
59 -1317.1797 1985.7559

> fit$CI[1,]
[1] -2191.498  2974.628
> range(BrachyBdhn1aD)
[1] 214.2 898.5
> obj<-contMap(BrachySpTreeD.tree,BrachyBdhn1aD,plot=FALSE)
> plot(obj,legend=0.7*max(nodeHeights(BrachySpTreeD.tree)),fsize=c(0.7,0.9))
