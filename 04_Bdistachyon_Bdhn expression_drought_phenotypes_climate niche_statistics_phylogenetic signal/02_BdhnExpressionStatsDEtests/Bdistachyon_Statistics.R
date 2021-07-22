##Script_by_Pilar_Catalan_26/04/2021 (Mod_MADR_JDC_30/4/21)
library(datasets)
library(stats)
library(base)
library(Rcpp)
library(lifecycle)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)
library(lattice)
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
tab_BrachyTPM_Bdhn<- read.csv(file="Bdistachyon_Bdhn_TPM_W&D.csv", h=T, sep=",", na.strings = "NA", stringsAsFactors=FALSE) 
tab_BrachyTPMall<- read.csv(file="Bdistachyon_Bdhn_Phenotypictraits_W&D.csv", h=T, sep=",", na.strings ="NA", stringsAsFactors=FALSE)
#Omit NAs
tab_all<- na.omit(tab_BrachyTPMall)  


#assign names to table columns
tab_all$ecotype<- tab_all[,1]
tab_all$Bdhn1a<- tab_all[,6]

#1. summarize main stats for each column
sink("Bdistachyon_TPM_W&D_MAIN_stats_all.csv")   
print(summary(tab_all))       
sink()                                

# 2. print histogram for each variable
hist(tab_BrachyTPMall$Bdhn1a)

# 2.1 boxplots of characters vs ecotypes 
png(filename="Boxplot_Bdhn1a.png", width = 1400, height = 788, units = "px")
boxplot(tab_BrachyTPM_Bdhn$Bdhn1a~ tab_BrachyTPM_Bdhn$ecotype,las=2, col=c("brown2","dodgerblue2"), main="Bdhn1a",
        xlab="", ylab="TPM")
dev.off()

# 3. Wilcoxon test for the whole table 
#remove unnnecessary columns
tab_all2<-subset(tab_all, select=-c(10))
tab_all3<- subset(tab_all2, select=-c(2:5))    
tab_all3$ecotype<-as.factor(tab_all3$ecotype)  #a factor  with the first column. It's use as groups
tab_all3.df<- data.frame(tab_all3)             
tab_all3.df$ecotype<-as.factor(tab_all3.df$ecotype)

levels(tab_all3.df$ecotype)

#3.1. Pairwise wilcox test for each variable in tab_all3 (done one by one)
#declaring variables
tab_all3.df$Bdhn1a<- as.numeric(as.character(tab_all3.df$Bdhn1a))

#Pairwise.wilcoxon test for an specific variable (preview step) against the whole set of variables
sink("pairwise.result.cn.csv")
pairwise.wilcox.test(tab_all3.df$Bdhn1a,tab_all3.df$ecotype,p.adjust.method = "BH")
sink()

#3.2. Wilcoxon pairwise test for each Bdhn by treatment and each Bdhn and phenotypic trait

contador=2  
while (contador <=16){
  versus<-colnames(subset( tab_all3,select=c(contador+1)))   
  Bdhn1aW<-unlist(subset( tab_all3, select= c(2)))         
  aux2<- unlist(subset( tab_all3, select= c(contador+1)))
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


# 3.3. Pairwise comparisons using Wilcoxon rank sum exact test 
# a data set was created for each ecotype
#Preparing dataset for each ecotype (slice for each set of row corresponding to an ecotype)
ABR2<- slice(tab_all2, 1:8)
ABR3<- slice(tab_all2,9:15)

#loop the preview dataset  for each ecotype
sink("ABR2_wilcox.csv")
contador=6
while (contador<=21){
  ABR2$contador<-ABR2[,contador]
  ABR2$treatment<-as.factor(ABR2$treatment)
  res <- wilcox.test( ABR2$contador~ ABR2$treatment , data = ABR2, exact = FALSE);res2<-res$p.value  
  print(contador)
  print(res)
  print(res2)
  contador=contador+1
}
sink()


# 4. kRUSKALL AND DUNNE
#create a variable joining the ecotype and treatment
tab_all2$ET <- paste(tab_all2$ecotype, tab_all2$treatment)
ET<-tab_all2$ET

# 4.1 Summarize by group
sink("summarize_taball2_ET.csv") 
contador=6
while(contador<=21){
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
delta13c_abs<-abs(delta13c)           
aux <- subset(tab_all2, select=-c(18))
tab_allABS<- cbind(aux, delta13c_abs)


contador=6
while(contador<=22){
  if (contador == 21) { contador=contador+1; next }
  Data1<- subset(tab_allABS, select=c(22,contador))
  tab_allABS$contador<-tab_allABS[,contador]
  hist<-histogram (~ tab_allABS$contador | tab_allABS$ecotype,
                   data=Data1,
                   layout=c(4,4), col="dodgerblue")      #  columns and rows of individual plot
  print(contador)
  print(hist)
  contador=contador+1
}



# 4.2 Kruskal-Wallis rank sum test
ET<- tab_allABS$ET
group = factor(ET, levels=unique(ET))

#4.2.1 Performing Kruskall and Dunn test in one step  

contador= 6
while(contador<=22){
  tab_allABS$contador<-tab_allABS[,contador]
  Krus<-kruskal.test(tab_allABS$contador, group) 
  Dunn<-posthoc.kruskal.dunn.test(tab_allABS$contador, group, data = tab_allABS, p.adjust="BH") #Dunn-test
  print(contador)
  print(Dunn)
  contador=contador+1
}



# 4.3. Pairwise de Wilcox, followed by Values tab, fulltab and multicomlletter for groups
sink("wilcox.pairwise_multlet.csv")
contador=6
while(contador<=22){
  if (contador == 21) { contador=contador+1; next }
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

#4.4. Analysis of lm, anova and sum the Variance Table (plot: plot_lm_an_tuk)
contador=6
while(contador<=22){
  if (contador == 21) { contador=contador+1; next }
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

# 4.5. Plotting Tukey 
dataD<-filter(tab_allABS, treatment=="Drought") #Drought filtered dataset
dataW<-filter(tab_allABS, treatment=="Watered") #Watered filtered dataset

#4.5.1 Drought
model=lm (dataD$Bdhn1a ~ dataD$ecotype)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dataD$ecotype', conf.level=0.95)
plot(TUKEY , las=1 , col="blue")

#Creating labels
generate_label_df <- function(TUKEY, variable){
  Tukey.levels <- TUKEY[[variable]][,4] # leave this value as 4. It is extracting the 4th column of the TUKEY object.
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'], stringsAsFactors = TRUE)
  Tukey.labels$ecotype=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$ecotype) , ]
  return(Tukey.labels)
}

LABELS <- generate_label_df(TUKEY , "dataD$ecotype")
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

b <- boxplot(dataD$Bdhn1a ~ dataD$ecotype , ylim=c(min(dataD$Bdhn1a) , 1.1*max(dataD$Bdhn1a)),
             las=2, col=my_colors[as.numeric(LABELS[,1])] , ylab="" , main="Bdhn1a")

over <- 0.2*max( b$stats[nrow(b$stats),] )

#Boxplots with Tukey letters

text( c(1:nlevels(factor(dataD$ecotype))) , b$stats[nrow(b$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )


#4.5.2 Watered
model=lm (dataW$Bdhn1a ~ dataW$ecotype)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dataW$ecotype', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")

#Creating labels
LABELS <- generate_label_df(TUKEY , "dataW$ecotype")
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
b <- boxplot(dataW$Bdhn1a ~ dataW$ecotype , ylim=c(min(dataW$Bdhn1a) , 1.1*max(dataW$Bdhn1a)) ,las=2, col=my_colors[as.numeric(LABELS[,1])] , ylab="%" , main="cn")

over <- 0.2*max( b$stats[nrow(b$stats),] )

#Boxplots with Tukey letters
text( c(1:nlevels(factor(dataD$ecotype))) , b$stats[nrow(b$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
