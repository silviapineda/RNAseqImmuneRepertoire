rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Analysis from the MIXCR data 
###           
### DESCRIP: Analysis with the data from the MIXCR output
###         
###
### Author: Silvia Pineda
### Date: March, 2018
############################################################################################
library("RColorBrewer")
library(plyr)
library(ggplot2)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/summaryMatrix.Rdata") ##Summary data from the MIXCR tool
load("Data/summaryCapture.Rdata") ##Summary data from the capture data in the QC process

##Put in same order both datasets
capture<-capture.df[order(rownames(capture.df)),]

##Total IG reads
summaryMatrix$IG_Reads<-summaryMatrix$IGH_Reads+summaryMatrix$IGK_Reads+summaryMatrix$IGL_Reads

###Analysis in Overall Ab expression
summaryMatrix$IG_expression<-summaryMatrix$IG_Reads/(summaryMatrix$IG_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$IGH_expression<-summaryMatrix$IGH_Reads/(summaryMatrix$IGH_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$IGK_expression<-summaryMatrix$IGK_Reads/(summaryMatrix$IGK_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$IGL_expression<-summaryMatrix$IGL_Reads/(summaryMatrix$IGL_Reads+as.numeric(capture$uniquelyMappedReads))

summaryMatrix$T_expression<-summaryMatrix$T_Reads/(summaryMatrix$T_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$TRA_expression<-summaryMatrix$TRA_Reads/(summaryMatrix$TRA_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$TRB_expression<-summaryMatrix$TRB_Reads/(summaryMatrix$TRB_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$TRD_expression<-summaryMatrix$TRD_Reads/(summaryMatrix$TRD_Reads+as.numeric(capture$uniquelyMappedReads))
summaryMatrix$TRG_expression<-summaryMatrix$TRG_Reads/(summaryMatrix$TRG_Reads+as.numeric(capture$uniquelyMappedReads))
clin<-factor(substr(rownames(capture),1,1))
clin<-revalue(clin, c("C"="CMR", "H"="AMR" , "N" = "STA"))
clin<-relevel(clin,ref="STA")

repertoireResults<-summaryMatrix
save(repertoireResults,clin,file="Data/RepertoireResults.Rdata")



COLOR=brewer.pal(3,"Set2")
tiff("Boxplot_IG_expression.tiff",res=300,w=2000,h=2000)
boxplot(summaryMatrix$IG_expression~clin,col=COLOR)
dev.off()

boxplot(summaryMatrix$IGH_expression~clin,col=COLOR)
boxplot(summaryMatrix$IGK_expression~clin,col=COLOR)
boxplot(summaryMatrix$IGL_expression~clin,col=COLOR)


tiff("Boxplot_T_expression.tiff",res=300,w=2000,h=2000)
boxplot(summaryMatrix$TRG_expression~clin,col=COLOR,main="T expression")
dev.off()
summary(glm(summaryMatrix$T_expression~clin)) #p(STA vs AMR)=0.03

tiff("Boxplot_TRs_expression.tiff",res=300,w=3000,h=3000)
par(mfrow=c(2,2))
boxplot(summaryMatrix$TRA_expression~clin,col=COLOR,main=c("TRA expression"))
boxplot(summaryMatrix$TRB_expression~clin,col=COLOR,main=c("TRB expression"))
boxplot(summaryMatrix$TRD_expression~clin,col=COLOR,main=c("TRD expression"))
boxplot(summaryMatrix$TRG_expression~clin,col=COLOR,main=c("TRG expression"))
dev.off()
summary(glm(summaryMatrix$TRA_expression~clin))
summary(glm(summaryMatrix$TRB_expression~clin))
summary(glm(summaryMatrix$TRD_expression~clin))
summary(glm(summaryMatrix$TRG_expression~clin))

tiff("Boxplot_alphabeta.tiff",res=300,w=2000,h=2000)
boxplot(summaryMatrix$Alpha_Beta~clin,col=COLOR,main="TRA+TRB/TRD+TRG ratio")
dev.off()
summary(glm(summaryMatrix$Alpha_Beta~clin)) #p(STA vs CMR)=0.01 p(STA vs AMR)=0.0001

plot_ratio<-cbind(data.frame(summaryMatrix$Alpha_Beta),clin)
colnames(plot_ratio)<-c("ratio","transplant_outcomes")
# Percent TRA & TRB of total T-reads
tiff(filename = "AlphaBeta.tiff", width = 5, height = 4, units = 'in', res = 300, compression = 'lzw')
ggplot(plot_ratio, aes(transplant_outcomes, ratio, fill = transplant_outcomes)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab(expression(paste(alpha,beta,"/",gamma,delta,"_ratio"))) + xlab("") + ylim(.85, 1)
dev.off()




