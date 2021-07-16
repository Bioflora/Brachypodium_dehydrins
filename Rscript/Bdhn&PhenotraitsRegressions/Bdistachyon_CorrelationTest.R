#Correlation Script Pilar Catalann, MADR-JDC(mayo2021) 

library(devtools)
library(readxl)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(psych)
library(PerformanceAnalytics)
library(ggplot2)
library(rcompanion)

#Import data
  # Bdistachyon general data
tab_all<-read.csv(file="Bdistachyon_Bdhn_PhenotypicTraitsW&D.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
tab_all$ecotype<-as.factor(tab_all$ecotype)
tab_all$Treatment<-as.factor(tab_all$Treatment)
tab_all$Climate<-as.factor(tab_all$Climate)
tab_all2<-na.omit(tab_all)
names(tab_all2)

  # Bdistachyon Watered and Drougth data split
tab_Bdis<-read.csv(file="Bdistachyon_Bdhns_PhenotypicTraits_WvsD.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE)
tab_Bdis2<-na.omit(tab_Bdis)
names(tab_Bdis2)


#1.histogram frequency
hist(log(tab_all2$rmr+1)) 

#1.1. Histogram by treatment (same for each variable)
hist(log(tab_Bdis2$rmrW+1)) 
hist(log(tab_Bdis2$rmrD+1))  


#2.linear regression model (lm) and estimates

##2.1 lm for each Bdhn  GENERAL (the loop must be run for each variable)
sink("lm_fit_Bdhn1a.csv")
contador2=5           #Must be change to run the loop for each variable
while(contador2<21){
  contador=6              
    while (contador <=22){
    aux_tab<-subset(tab_all2, select=-c(10)) #delete treatment column.
    versus<-colnames(subset( aux_tab,select=c(contador+1)))  
    versus1<-colnames(subset(aux_tab, select=c(contador2+1)))
    aux1<-unlist(subset( aux_tab, select= c(contador2+1)))
    aux2<- unlist(subset( aux_tab, select= c(contador+1)))
    plot(aux1~aux2, main=paste(versus), ylab=versus1)
    fit<-summary(fit<-lm(aux1~aux2));abline(fit,col="darkgreen",lty=1, lwd=2)
    print(fit)                  
    contador = contador+1      
    }
  contador2=contador2+1
  
}
sink()


# 2.2 linear model for each Bdhn by treatment 
  #plot Bdhn1aW vs Bdhn1aD linear regression model (lm) and estimates

#Watered loop
sink("fit_Bdhn1aW.csv")
contador=2    
while (contador <=13){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   
  Bdhn7W<-unlist(subset( tab_Bdis2, select= c(5)))           #Change name
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  plot(Bdhn7W~aux2, main=paste(versus))
  fit<-summary(fit<-lm(Bdhn7W~aux2));abline(fit,col="blue",lty=1, lwd=2)
  print(plot)
  print(fit)                  
  contador = contador+1       
}
sink()


#Drought loop 
sink("fit_Bdhn1aD.csv")
contador=18                #must be incremented and run the loop 4 times.
while (contador <=33){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   
  Bdhn1aD<-unlist(subset( tab_Bdis2, select= c(20)))          #Change name
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  plot(Bdhn1aD~aux2, main=paste(versus))
  fit<-summary(fit<-lm(Bdhn1aD~aux2));abline(fit,col="orange",lty=1, lwd=2)
  print(fit)                  
  contador = contador+1      
}
sink()

#3. Wilcoxon pairwise test for each Bdhn by treatment and each Bdhn and phenotypic trait

contador=2  
while (contador <=32){
  versus<-colnames(subset( tab_Bdis2,select=c(contador+1)))   
  Bdhn1aW<-unlist(subset( tab_Bdis2, select= c(2)))         
  aux2<- unlist(subset( tab_Bdis2, select= c(contador+1)))
  Wilcox<-wilcox.test( aux2,Bdhn1aW, paired=TRUE)
  png(paste("Wilcox_Bdhn1a", versus, ".png", sep = ""), width=955, height= 450)  
  plot(aux2, Bdhn1aW,
       pch = 16,
       main=versus)
  abline(0,1, col="blue", lwd=2)
  dev.off()
  print(Wilcox)
  contador = contador+1 
  }


# Data for Correlation test
tab_Bdhns <- subset(tab_Bdis2, select=c(2:5, 17:21 ))

tab_Bdhns$Bdhn1aW<-as.numeric(tab_Bdhns$Bdhn1aW)
tab_Bdhns$Bdhn2W<-as.numeric(tab_Bdhns$Bdhn2W)
tab_Bdhns$Bdhn3W<-as.factor(tab_Bdhns$Bdhn3W)
tab_Bdhns$Bdhn7W<-as.factor(tab_Bdhns$Bdhn7W)
tab_Bdhns$Bdhn1aD<-as.factor(tab_Bdhns$Bdhn1aD)
tab_Bdhns$Bdhn2D<-as.factor(tab_Bdhns$Bdhn2D)
tab_Bdhns$Bdhn3D<-as.factor(tab_Bdhns$Bdhn3D)
tab_Bdhns$Bdhn7D<-as.factor(tab_Bdhns$Bdhn7D)



# 4 Correlacion tests (Pearson, Spearman & Kendall)

  #4.1 Function to extract correlation coefficient and p-values
corrFuncP <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2], method="pearson", exact=FALSE)
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
             stringsAsFactors=FALSE)
}

#The loop to run Pearson correlation test in all of the variables.
  
  #Pearson Test
sink("Correlation_testPearson.csv", type="output")
contador=2
while(contador<=32){
vars2 = data.frame(v1=names(tab_Bdis2)[contador], v2=names(tab_Bdis2)[-1]) # Apply corrFunc to all rows of vars

P1 = do.call(rbind, mapply(corrFuncP, vars2[,(1)], vars2[,2], MoreArgs=list(data=tab_Bdis2), 
                                 SIMPLIFY=FALSE))
print(P1)
contador=contador+1
}
sink()

 

# 5 multivariate correlation (https://rcompanion.org/handbook/I_10.html)
  #5.1 Correlation between all of the variables

tab_all$ecotype = factor(tab_all$ecotype,
                         levels=unique(tab_all$ecotype))
summary(tab_all)

pairs(data=tab_all,
      ~ Bdhn1a + Bdhn2 + Bdhn3 + Bdhn7 + leaf_rwc + leaf_wc+ lma + pro + abvrgd + 
        blwgrd +ttlmass + rmr + delta13c + leafc + leafn +cn ) 


Data.num = tab_all[c("Bdhn1a", "Bdhn2", "Bdhn3", "Bdhn7", "leaf_rwc","leaf_wc", "lma", "pro", "abvrgd", "blwgrd",
                  "ttlmass", "rmr",  "delta13c", "leafc", "leafn", "cn")]
print(corr.test(Data.num,
          use    = "pairwise",
          method = "pearson",
          adjust = "none"),
      short=FALSE)  


  #5.2 Correlation between all of the variables for each Bdhn

tabB4<-subset(tab_all, select= -c(2:8, 10 )) #Make a subset for each Bdhn

tabB4$ecotype = factor(tabB4$ecotype, levels=unique(tabB4$ecotype))
summary(tabB4)

pairs(data=tabB4,
      ~ Bdhn7 + leaf_rwc + leaf_wc+ lma + pro + abvrgd + 
        blwgrd +ttlmass + rmr + delta13c + leafc + leafn +cn ) 

Data.num2<-subset(tabB4, select= -c(1))
corr.test(Data.num2,
          use    = "pairwise",
          method = "pearson",
          adjust = "none")

