#Script Pilar Catalán (), MADR-JDC(mayo2021) 

#correlation scrips
install.packages("ggpubr")
library(devtools)
library(readxl)
library(ggplot2)
library(ggsignif)
library(ggpubr)


tab_all<-read.csv(file="BrachyBdhnPhenoall.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
tab_all1<- na.omit(tab_all)

tab_all$ecotype<-as.factor(tab_all$ecotype)
tab_all$Treatment<-as.factor(tab_all$Treatment)
tab_all$Climate<-as.factor(tab_all$Climate)
names(tab_all)

tab_Bdis<-read.csv(file="BrachyBdhnPhenoallWvsDmod.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
names(tab_Bdis)


######################################
  #Sergios's
#sum of columns 8 to 33 into a new column 'sum'
#Shannon diversity index for columns 8 to 33 into a new column H_vocs
#no.of VOC columns with data no zero('species number') of columns 8 to 33 into a new column S
#omit all NA in all columns (-13?)

tab_all$sum<-rowSums(tab_all[,8:33])
tab_all$H_vocs<-diversity(tab_all[,8:33])
tab_all$S <- specnumber(tab_all[,8:33])

##################################################

#omit all NA in all columns
tab_all2<-na.omit(tab_all)
tab_Bdis2<-na.omit(tab_Bdis)

#1.histogram of Bdhn1a frequency
hist(log(tab_all2$rmr+1)) 

hist(log(tab_Bdis2$rmrW+1))  #idem para tab_bdis2
hist(log(tab_Bdis2$rmrD+1))  

#2. plot Bdhn1a vs Bdhn2 linear regression model (lm) and estimates (ver while)

plot(Bdhn1a~rmr, tab_all2)

##2.1 While con contador para lm de cada Bdhn by MADR_JDC   GENERAL

contador2=5
while(contador2<9){
  contador=5
    while (contador <=9){
    aux_tab<-subset(tab_all2, select=-c(10)) #quita la columna 10 (treatment)
    versus<-colnames(subset( aux_tab,select=c(contador+1)))  #Crea una variable para el titulo segun va corriendo el contador
    versus1<-colnames(subset(aux_tab, select=c(contador2+1)))
    aux1<-unlist(subset( aux_tab, select= c(contador2+1)))           #Para la columna 6 = Bdhn1a (correrx4)
    aux2<- unlist(subset( aux_tab, select= c(contador+1)))
    plot(aux1~aux2, main=paste(versus))
    fit<-summary(fit<-lm(aux1~aux2));abline(fit,col="blue",lty=1, lwd=2)
    print(fit)                  #Permite que se imprima en pantalla el fit
    contador = contador+1       #cambiar de columnas en el loop
    }
  contador2=contador2+1
  
}

######################
sink("lm_fit_Bdhn7.csv")
contador2=8
while(contador2<=10){
  contador=8
  while (contador <=22){
    aux_tab<-subset(tab_all2, select=-c(10)) #quita la columna 10 (treatment)
    versus<-colnames(subset( aux_tab,select=c(contador+1)))  #Crea una variable para el titulo segun va corriendo el contador
    versus1<-colnames(subset(aux_tab, select=c(contador2+1)))
    auxAA<-unlist(subset( aux_tab, select= c(contador2+1)))           #Para la columna 6 = Bdhn1a (correrx4)
    aux2<- unlist(subset( aux_tab, select= c(contador+1)))
    plot(auxAA~aux2, main=paste(versus), xlab=versus1, ylab=versus)
    fit<-summary(fit<-lm(auxAA~aux2));abline(fit,col="chartreuse4",lty=1, lwd=2)
    print(versus)
    print(fit)                  #Permite que se imprima en pantalla el fit
    contador = contador+1       #cambiar de columnas en el loop
  }
  contador2=contador2+1
}
sink()
##############################################################

sink("plots_pheno.csv")
contador=11
while (contador <=22){
 aux_tab1<-subset(tab_all2, select=-c(10))                 #quita la columna 10 (treatment)
 versus<-colnames(subset( aux_tab,select=c(contador+1)))  #Crea una variable para el titulo segun va corriendo el contador
 leafrwc<-unlist(subset( aux_tab, select= c(11)))           #Para la columna 6 = Bdhn1a (correrx4)
 aux2<- unlist(subset( aux_tab, select= c(contador+1)))
 plot(leafrwc~aux2, main=paste(versus))
 fit<-summary(fit<-lm(leafrwc~aux2));abline(fit,col="green",lty=1, lwd=2)
 print(fit)  #Permite que se imprima en pantalla el fit
 print(plot)
 contador = contador+1       #cambiar de columnas en el loop
 }

sink()




#################################################################################################################
# 2.2 linear model for each Bdhn by treatment 
#plot Bdhn1aW vs Bdhn1aD linear regression model (lm) and estimates

sink("fit_conRMRWBdhn7.csv")
contador=2    #Para W
while (contador <=13){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   #Crea una variable para el titulo segun va corriendo el contador
  Bdhn7W<-unlist(subset( tab_Bdis2, select= c(5)))           #Para la columna 2 = Bdhn1aW (correrx4)
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  plot(Bdhn7W~aux2, main=paste(versus))
  fit<-summary(fit<-lm(Bdhn7W~aux2));abline(fit,col="blue",lty=1, lwd=2)
  print(plot)
  print(fit)                  #Permite que se imprima en pantalla el fit
  contador = contador+1       #cambiar de columnas en el loop
}
sink()

####################


Bdhn1aD<-unlist(subset( tab_Bdis2, select= c(18)))           #Para la columna 2 = Bdhn1aW (correrx4)
ProD<- unlist(subset( tab_Bdis2, select= c(25)))
plot(Bdhn1aD~ProD, main=paste("proD"))
fit<-summary(fit<-lm(Bdhn1aD~ProD));abline(fit,col="blue",lty=1, lwd=2)
print(plot)
print(fit)                  #Permite que se imprima en pantalla el fit





#plot Bdhn1aD vs Bdhn1aD linear regression model (lm) and estimates

contador=17  #Para D
while (contador <=33){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   #Crea una variable para el titulo segun va corriendo el contador
  Bdhn1aD<-unlist(subset( tab_Bdis2, select= c(20)))           #Para la columna var
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  plot(Bdhn7D~aux2, main=paste(versus))
  fit<-summary(fit<-lm(Bdhn1aD~aux2));abline(fit,col="orange",lty=1, lwd=2)
  print(fit)                  #Permite que se imprima en pantalla el fit
  contador = contador+1      #cambiar de columnas en el loop
}

####################################################################################################################

# 3. wilcoxon pairwise test  #para ccada Bdhn en cada tratamiento y cada Bdhn x pheno

contador=2  #Test wilcox
while (contador <=32){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   #Crea una variable para el titulo segun va corriendo el contador
  Bdhn1aW<-unlist(subset( tab_Bdis2, select= c(2)))           #Para la columna 2 = Bdhn1aW (correrx4)
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  Wilcox<-wilcox.test( aux2,Bdhn1aW, paired=TRUE)
  #png(paste("Wilcox_Bdhn1a", versus, ".png", sep = ""), width=955, height= 450)  ##Para imprimir en jpg en la carpeta directorio trabajo
  plot(aux2, Bdhn1aW,
       pch = 16,
       main=versus)
  abline(0,1, col="blue", lwd=2)
  #dev.off()
  print(Wilcox)
  contador = contador+1 #cambiar de columnas en el loop
  }

Bdhn1aW<-unlist(subset( tab_Bdis2, select= c(2)))           #Para la columna 2 = Bdhn1aW (correrx4)
leaf_wc<- unlist(subset( tab_Bdis2, select= c(7)))
Wilcox<-wilcox.test( leaf_wc,Bdhn1aw, paired=TRUE)

#Plot of paired samples from a Wilcoxon signed-rank test.
#Circles above and to the left of the blue one-to-one line indicate observations with a higher value


####Test de correlacion ########

tab_Bdhns <- subset(tab_Bdis2, select=c(2:5, 17:20 ))
tab_Bdhns
tab_Bdhns$Bdhn1aW<-as.numeric(tab_Bdhns$Bdhn1aW)
tab_Bdhns$Bdhn1aW
tab_Bdhns$Bdhn2W<-as.numeric(tab_Bdhns$Bdhn2W)
tab_Bdhns$Bdhn2W
tab_Bdhns$Bdhn3W<-as.factor(tab_Bdhns$Bdhn3W)
tab_Bdhns$Bdhn7W<-as.factor(tab_Bdhns$Bdhn7W)
tab_Bdhns$Bdhn1aD<-as.factor(tab_Bdhns$Bdhn1aD)
tab_Bdhns$Bdhn2D<-as.factor(tab_Bdhns$Bdhn2D)
tab_Bdhns$Bdhn3D<-as.factor(tab_Bdhns$Bdhn3D)
tab_Bdhns$Bdhn7D<-as.factor(tab_Bdhns$Bdhn7D)


###########################
#test de correlacion
	#este bloque permite formatear los resultados en una tabla
# Function to extract correlation coefficient and p-values

#4.1 Correlacion test (Pearson, Spearman & Kendall)
corrFuncP <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2], method="pearson", exact=FALSE)
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
             stringsAsFactors=FALSE)
}

corrFuncS <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2], method="spearman", exact=FALSE)
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
             stringsAsFactors=FALSE)
}
corrFuncK <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2], method="kendall", exact=FALSE)
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
             stringsAsFactors=FALSE)
}


    # Pairs of variables for which we want correlations
#vars = data.frame(v1=names(tab_Bdhns)[1], v2=names(tab_Bdhns)[-1])

    # Apply corrFunc to all rows of vars
#corrs_P = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=tab_Bdhns), 
#                              SIMPLIFY=FALSE))
sink("Correlation_testKendall.csv", type="output")
contador=2
while(contador<=32){
vars2 = data.frame(v1=names(tab_Bdis2)[contador], v2=names(tab_Bdis2)[-1])

#P1 = do.call(rbind, mapply(corrFuncP, vars2[,(1)], vars2[,2], MoreArgs=list(data=tab_Bdis2), 
#                                 SIMPLIFY=FALSE))

K1 = do.call(rbind, mapply(corrFuncK, vars2[,1], vars2[,2], MoreArgs=list(data=tab_Bdis2), 
                                 SIMPLIFY=FALSE))

print(K1)
contador=contador+1
}
sink()


### ¿?
ggscatter(tab_Bdhns, x = "Bdhn1aW", y = "Bdhn1aD", 
          add = "reg.line", cor.coef = TRUE,
          conf.int = TRUE, cor.method = "kendall",
          xlab = "Bdhn1aW", ylab = "Bdhn1aD", main= "Spearman", color="blue")



####################################################################################################################

# multivariate correlation (https://rcompanion.org/handbook/I_10.html)

if(!require(psych)){install.packages("psych")}
if(!require(PerformanceAnalytics)){install.packages("PerformanceAnalytics")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(rcompanion)){install.packages("rcompanion")}
library(psych)
library(PerformanceAnalytics)
library(ggplot2)
library(rcompanion)

tab_all1$ecotype = factor(tab_all1$ecotype,
                         levels=unique(tab_all1$ecotype))
summary(tab_all1)

pairs(data=tab_all1,
      ~ Bdhn1a + Bdhn2 + Bdhn3 + Bdhn7 + leaf_rwc + leaf_wc+ lma + pro + abvrgd + 
        blwgrd +ttlmass + rmr + delta13c + leafc + leafn +cn ) #enumera las variables a comparar


Data.num = tab_all1[c("Bdhn1a", "Bdhn2", "Bdhn3", "Bdhn7", "leaf_rwc","leaf_wc", "lma", "pro", "abvrgd", "blwgrd",
                  "ttlmass", "rmr",  "delta13c", "leafc", "leafn", "cn")]
print(corr.test(Data.num,
          use    = "pairwise",
          method = "pearson",
          adjust = "none"),
      short=FALSE)  #si uso print (corr(), short=FALSE) se imprimen los intervalos de confianza


chart.Correlation(Data.num,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)

##########

tabB4<-subset(tab_all1, select= -c(2:8, 10 ))

tabB4$ecotype = factor(tabB4$ecotype, levels=unique(tabB4$ecotype))
summary(tabB4)

pairs(data=tabB4,
      ~ Bdhn7 + leaf_rwc + leaf_wc+ lma + pro + abvrgd + 
        blwgrd +ttlmass + rmr + delta13c + leafc + leafn +cn ) #enumera las variables a comparar

#Data.num2 = tabB1a[c("Bdhn1a", "leaf_rwc","leaf_wc", "lma", "pro", "abvgrd", "blwgrd",
 #                     "ttlmass", "rmr",  "delta13c", "leafc", "leafn", "cn")]  #se supone que hace un subset pero ahora no

Data.num2<-subset(tabB4, select= -c(1))
corr.test(Data.num2,
          use    = "pairwise",
          method = "pearson",
          adjust = "none")

chart.Correlation(Data.num2,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)

