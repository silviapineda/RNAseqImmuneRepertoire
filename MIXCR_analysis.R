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




