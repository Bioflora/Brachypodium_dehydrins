#Script_by_Pilar_Catalan_26/04/2021 (Mod_MADR_JDC_30/4/21)
update.packages(checkBuilt = TRUE)
install.packages("Rcpp")
install.packages("plotly")
install.packages("lifecycle")
install.packages("ellipsis")
install.packages("dplyr")
install.packages("stats")
install.packages("base")
install.packages("DBI")
install.packages("PMCMRplus")
install.packages("lapply")
library(datasets)
library(stats)
library(base)
library(Rcpp)
library(lifecycle)
library(tidyverse)
library(dplyr)
library("ggpubr")#para los wilcox whiskers plot
library(lattice)
library("lapply")
library(FSA)
library(caret)
library(rcompanion)
library(multcompView)
library(car)
library(agricolae)
library(PMCMR)
library(PMCMRplus)
library(ggplot2)
library(ggsignif)




#import data as table
tab_BrachyTPM_Bdhn<- read.csv(file="BrachyTPM_W&D_stats2.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE) #solo las Bdhns
tab_BrachyTPMall<- read.csv(file="BrachyBdhnPhenoW&Dallmod2prueba.csv", h=T, sep=",", na.strings ="NA", stringsAsFactors=FALSE)


tab_BrachyTPMall_1<- na.omit(tab_BrachyTPMall)  #Tabla completa pero sin NANs

#inspect table columns
ecotype <- tab_BrachyTPMall_1 [,1]
Bdhn1a <- tab_BrachyTPMall_1 [,6]
leaf_rwc<-tab_BrachyTPMall_1[, 11]

#assign names to table columns
tab_BrachyTPMall_1$ecotype<- tab_BrachyTPMall_1[,1]
tab_BrachyTPMall_1$Bdhn1a<- tab_BrachyTPMall_1[,6]



#visualize the data, plots of all pairwise columns' values and save it at png
png( filename = "tab_BrachyTPMall.png", width = 491, height = 739)
plot(tab_BrachyTPMall_1)
dev.off()
#plot(tab_BrachyTPMall_1$Bdhn1a, tab_BrachyTPMall_1$leaf_rwc) #plot a species pairwise comparison of two variables


# 1. summarize main stats for each column

sink("BrachyTPM_W&D_MAIN_stats_all.csv")   
print(summary(tab_BrachyTPMall_1))       
sink()                                


# 2.1. print histogram for each variable
hist(tab_BrachyTPMall$Bdhn7)

# 2.2. boxplots of characters vs ecotypes (las=2 to add names/numerical scales in x and y axes; col= for specific colors)
png(filename="Boxplot_Bdhn1.png", width = 1400, height = 788, units = "px")
boxplot(tab_BrachyTPM_Bdhn$Bdhn1a~ tab_BrachyTPM_Bdhn$ecotype,las=2, col=c("brown2","dodgerblue2"), main="Bdhn1a",
        xlab="", ylab="TPM")
dev.off()


# 2.3.  boxplots with rainbow colors (for multiple samples)
c1 <- rainbow(10)
c2 <- rainbow(10, alpha=0.2)
c3 <- rainbow(10, v=0.7)
boxplot(tab_BrachyTPMall$Bdhn3 ~ tab_BrachyTPMall$ecotype,las=2,col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=2)



# 3. Wilcoxon test for the whole table MD

tab_all<- tab_BrachyTPMall_1

#remove unnnecessary columns
tab_all2<-subset(tab_all, select=-c(10))
tab_all3<- subset(tab_all2, select=-c(2:5))    #next step only
tab_all3$ecotype<-as.factor(tab_all3$ecotype)  #a factor  with the first column. It's use as groups
tab_all3.df<- data.frame(tab_all3)             
tab_all3.df$ecotype<-as.factor(tab_all3.df$ecotype)

levels(tab_all3.df$ecotype)

#3.1. Pairwise wilcox test for each variable in tab_all3 (done one by one)
   #declaring variables
tab_all3.df$Bdhn1a<- as.numeric(as.character(tab_all3.df$Bdhn1a))
tab_all3.df$Bdhn2<- as.numeric(as.character(tab_all3.df$Bdhn2))


sink("pairwise.result.cn.csv")
pairwise.wilcox.test(tab_all3.df$leafn,tab_all3.df$ecotype,p.adjust.method = "BH")
sink()


# 3.2. Pairwise comparisons using Wilcoxon rank sum exact test 
       # a data set was created for each ecotype
      
#extraccion de los dataset de cada ecotipo mediante slice(#row)
ABR2<- slice(tab_all2, 1:8)
ABR3<- slice(tab_all2,9:15)

#loop for each ecotype
sink("Bd21_3wilcox.csv")
contador=6
while (contador<=21){
  Bd21_3$contador<-Bd21_3[,contador]
  Bd21_3$treatment<-as.factor(Bd21_3$treatment)
  res <- wilcox.test( Bd21_3$contador~ Bd21_3$treatment , data = Bd21_3, exact = FALSE);res2<-res$p.value  
  print(contador)
  print(res)
  print(res2)
  contador=contador+1
}
sink()


# 4. kRUSKALL AND DUNNE
  #create a variable joining the ecoyte and treatment
tab_all2$ET <- paste(tab_all2$ecotype, tab_all2$treatment)
tab_all2$ET


# 4.1 Summarize by group
contador=6
sink("summarize_taball2_ET_12may.csv")
while(contador<=6){
  Data1<- subset(tab_all2, select=c(22,contador))
    nombre<-unlist(subset( tab_all2,select=c(contador)))
  Data2 = mutate(Data1,
               group = factor(ET, levels=unique(ET))) 
  sum<-Summarize(nombre~group,
         data=Data2)
  print(contador)
  print(sum)
contador=contador+1
}
sink()


#4.1.1Plotting
  #delta13c is negative so it must be positive to be plotted
tab_all2$ET <- paste(tab_all2$ecotype, tab_all2$treatment)
delta13c= as.data.frame(tab_all2 [18, drop=FALSE])
delta13c_abs<-abs(delta13c)             #absolute values for this variable
aux <- subset(tab_all2, select=-c(18))
tab_allABS<- cbind(aux, delta13c_abs)

png( filename = paste("tab_BrachyTPM_var", contador, ".png", sep=""), width = 491, height = 739) 
plot(tab_BrachyTPMall_1)
dev.off()



contador=6
while(contador<=6){
  Data1<- subset(tab_allABS, select=c(22,contador))
  tab_allABS$contador<-tab_allABS[,contador]
  hist<-histogram (~ tab_allABS$contador | ecotype,
          data=Data1,
          layout=c(4,4), col="dodgerblue")      #  columns and rows of individual plot
  print(contador)
  print(hist)
  contador=contador+1
}


#4.2. Performing Kruskall and Dunn test in one step
# 4.3 Kruskal-Wallis rank sum test
tab_allABS$Bdhn2<-tab_allABS[,7]
aux<-tab_allABS$Bdhn3<-tab_allABS[,8]
tab_allABS$ET<-tab_allABS[,22]

ET<- tab_allABS$ET
group = factor(ET, levels=unique(ET))
group
K<-kruskal.test(tab_allABS$Bdhn1a, group) 
K
D<-posthoc.kruskal.dunn.test(tab_allABS$Bdhn2 ~ group, data = tab_allABS, p.adjust="BH") #Dunn-test

sink("krus_Dun.csv")
contador= 6
while(contador<=6){
  tab_allABS$contador<-tab_allABS[,contador]
  Krus<-kruskal.test(tab_allABS$contador, group) 
  Dunn<-posthoc.kruskal.dunn.test(tab_allABS$contador, group, data = tab_allABS, p.adjust="BH") #Dunn-test
  print(contador)
  print(Dunn)
  contador=contador+1
}

sink()


#Dunn test for multiple comparisons --> hecho en el comando anterior
#If the Kruskal?Wallis test is significant, a post-hoc analysis can be performed 
#to determine which levels of the independent variable differ from each other level. 


# pairwise.wilcoxon (Mann?Whitney U-tests.)
#Another post-hoc approach is to use pairwise Mann?Whitney U-tests.  
#To prevent the inflation of type I error rates, adjustments to the p-values can 
#be made using the p.adjust.method option to control the familywise error rate or
#to control the false discovery rate. See ?p.adjust for details.
#If there are several values to compare, it can be beneficial to have R convert
#this table to a compact letter display for you.  The multcompLetters function in 
#the multcompView package can do this, but first the table of p-values must be converted to a full table.


# 4.4. Pairwise de Wilcox, followed by Values tab, fulltab and multicomlletter for groups
sink("wilcox.pairwise_multlet.csv")
contador=6
while(contador<=6){
  tab_allABS$contador<-tab_allABS[,contador]  
  PT3 = pairwise.wilcox.test(tab_allABS$contador,
                             group,
                             p.adjust.method="BH")
  PT4 = PT3$p.value
  PT5 = fullPTable(PT4)
  multcompletters<- multcompLetters(PT4,
                                    compare="<",
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE)
 
  
  print(contador)
  print(PT3)
  print(PT4)
  print(PT5)
  print(multcompletters)
  contador=contador+1
}
sink()

#4.5. Analysis of lm, anova and sum the Variance Table (plot: plot_lm_an_tuk)

contador=6
while(contador<=6){
  tab_allABS$contador<-tab_allABS[,contador]
  model = lm(tab_allABS$contador ~ group,
           data=tab_allABS)
    hist(residuals(model),
       col="dodgerblue")
    a<-Anova(model, type="II")
    summary(model)
    plot(fitted(model),
       residuals(model))
  TUKEY.test<-(HSD.test(model, "group"))
  print(contador)
  print(Anova(model, type="II"))
  print(anova(model))
  print(summary(model))
  print(TUKEY.test)
  contador=contador+1
}

#########################################################################################################################

# 4.6. Plotting Tukey 
#Para drought

library(tidyverse)
data<-filter(tab_allABS, treatment=="Drought") #filtrado para Drougth
write.table(data, file="dataT.csv", append=FALSE, quote=TRUE, sep=",") #Extraido para quitar los guiones
data<- read.csv(file="dataTuk.csv", h=T, sep=";", na.strings = "NA", stringsAsFactors=TRUE) #introduce nueva
data
#PAra watered
dataW<-filter(tab_allABS, treatment=="Watered") #filtrado para Drougth
dataW
write.table(dataW, file="dataW.csv", append=FALSE, quote=TRUE, sep=",") #Extraido para quitar los guiones
dataW1<- read.csv(file="dataW1.csv", h=T, sep=";", na.strings = "NA", stringsAsFactors=TRUE) #introduce nueva
dataW1
###################################################################################################################

#Drought
model=lm (data$cn ~ data$ecotype)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'data$ecotype', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")


#library(multcompView) # load this library
generate_label_df <- function(TUKEY, variable){
  
  
  Tukey.levels <- TUKEY[[variable]][,4] # leave this value as 4. It is extracting the 4th column of the TUKEY object.
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'], stringsAsFactors = TRUE)
  
  
  Tukey.labels$ecotype=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$ecotype) , ]
  return(Tukey.labels)
}
LABELS <- generate_label_df(TUKEY , "data$ecotype")
my_colors <- c( 
  rgb(178, 34, 34, maxColorValue = 255), 
  rgb(255, 0, 10,    maxColorValue = 255), 
  rgb(255, 140, 0,  maxColorValue = 255), 
  rgb(253, 142, 36  ,maxColorValue = 255), 
  rgb(253, 206, 36 , maxColorValue = 255), 
  rgb(202, 18, 18  ,maxColorValue = 255), 
  rgb(202, 18, 18 , maxColorValue = 255), 
  rgb(202, 18, 18  ,maxColorValue = 255), 
  rgb(202, 18, 18  ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18  ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18  ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18  ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18  ,maxColorValue = 255),
  rgb(202, 18, 18   ,maxColorValue = 255),
  rgb(202, 18, 18  ,maxColorValue = 255)
  )
b <- boxplot(data$cn ~ data$ecotype , ylim=c(min(data$cn) , 1.1*max(data$cn)) ,las=2, col=my_colors[as.numeric(LABELS[,1])] , ylab="" , main="cn")

over <- 0.02*max( b$stats[nrow(b$stats),] )

text( c(1:nlevels(data$ecotype)) , b$stats[nrow(b$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )



#Watered
model=lm (dataW1$cn ~ dataW1$ecotype)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dataW1$ecotype', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")


#library(multcompView) # load this library
generate_label_df <- function(TUKEY, variable){
  
  
  Tukey.levels <- TUKEY[[variable]][,4] # leave this value as 4. It is extracting the 4th column of the TUKEY object.
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'], stringsAsFactors = TRUE)
  
  
  Tukey.labels$ecotype=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$ecotype) , ]
  return(Tukey.labels)
}
LABELS <- generate_label_df(TUKEY , "dataW1$ecotype")
my_colors <- c( 
  rgb(15, 42, 159, maxColorValue = 255), 
  rgb(30, 144, 255, maxColorValue = 255), 
  rgb(0, 191, 255, maxColorValue = 255), 
  rgb(70, 130, 180, maxColorValue = 255), 
  rgb(15, 42, 159, maxColorValue = 255), 
  rgb(15, 42, 159, maxColorValue = 255), 
  rgb(15, 42, 159, maxColorValue = 255), 
  rgb(15, 42, 159, maxColorValue = 255), 
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255),
  rgb(15, 42, 159, maxColorValue = 255)
  
)
b <- boxplot(dataW1$cn ~ dataW1$ecotype , ylim=c(min(dataW1$cn) , 1.1*max(dataW1$cn)) ,las=2, col=my_colors[as.numeric(LABELS[,1])] , ylab="%" , main="cn")

over <- 0.02*max( b$stats[nrow(b$stats),] )

text( c(1:nlevels(dataW1$ecotype)) , b$stats[nrow(b$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )

