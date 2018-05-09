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

load("Data/repertoireResults.Rdata")

##### Analysis of the different reads from MIXCR to find association with clinical outcome
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
summary(glm(summaryMatrix$TRA_expression~clin)) #(STA vs AMR p=0.002)
summary(glm(summaryMatrix$TRB_expression~clin)) #(STA vs AMR p=0.03)
summary(glm(summaryMatrix$TRD_expression~clin)) #(STA vs CMR p=0.00008 STA vs AMR p=0.0002)
summary(glm(summaryMatrix$TRG_expression~clin)) #(STA vs CMR p=0.06 STA vs AMR p=0.01)

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


##### Clinical Data Analysis #######
clin_data<-read.csv("Data/ClinicalData.csv")

id_clin<-match(rownames(summaryMatrix),clin_data$Individual_id)
clin_data<-clin_data[id_clin,]
clin_data$clin<-clin

summary(glm(summaryMatrix$TRG_expression~clin_data$DSA_HLA.class..0.neg..1.class.I..2.class.II..3.class.I.and.II.)) #
boxplot(summaryMatrix$Alpha_Beta~clin_data$DSA_HLA.class..0.neg..1.class.I..2.class.II..3.class.I.and.II.,col=COLOR)



#####################
#### GTEX data ######
####################
load("Data/SummaryMatrixReadsFromMIXCR_GTEX.Rdata")

ratio<-c(summaryMatrix$Alpha_Beta_ratio_expression,summaryMatrix_Everything$AlphaBeta_Percentage)
clin2<-factor(c(as.character(clin),rep("GTEx",length(summaryMatrix_Everything$AlphaBeta_Percentage))))
summary(glm(ratio~relevel(clin2,ref="GTEx"))) #
boxplot(ratio~clin2,col=COLOR)


gtex_id<-read.csv("Data/GTEX/SraRunTable_blood_1691.csv")
individual_id<-gtex_id[match(rownames(summaryMatrix_Everything),gtex_id$Run_s),"submitted_subject_id_s"]
gtex_phenotype_data<-read.csv("Data/GTEX/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.csv")

###clasify diseases by renal, others and unknown
gtex_phenotype_data$causeOfDeath<-ifelse(gtex_phenotype_data$DTHCOD=="acute renal failure" | 
gtex_phenotype_data$DTHCOD=="acute renal failure secondary to polycystic kidney disease" | gtex_phenotype_data$DTHCOD=="end-stage renal disease"
| gtex_phenotype_data$DTHCOD=="ESRD" | gtex_phenotype_data$DTHCOD=="esrd (end stage renal disease)" | gtex_phenotype_data$DTHCOD=="failure, renal" 
| gtex_phenotype_data$DTHCOD=="Kidney diseases" | gtex_phenotype_data$DTHCOD=="kidney failure" | gtex_phenotype_data$DTHCOD=="renal failure" | gtex_phenotype_data$DTHCOD=="Renal failure","kidney-related",
ifelse(gtex_phenotype_data$DTHCOD=="death - cause unknown" | gtex_phenotype_data$DTHCOD=="unknown" | gtex_phenotype_data$DTHCOD=="unknown cause of death" 
       | gtex_phenotype_data$DTHCOD=="unknown, do not have a copy of the death certificate or an ME report" 
       | gtex_phenotype_data$DTHCOD=="Unknown, do not have a copy of the death certificate or an ME report", "Unknown","Other-causes"))

COD<-gtex_phenotype_data[match(individual_id,gtex_phenotype_data$SUBJID),"causeOfDeath"]
clin3<-factor(c(as.character(clin),as.character(COD)))

summary(glm(ratio~relevel(clin3,ref="STA"))) #
boxplot(ratio~clin3,col=COLOR)
