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
library(gridExtra)
library(grid)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/repertoireResults.Rdata")

##### Analysis of the different reads from MIXCR to find association with clinical outcome

###Summary of the IG's
mean(summaryMatrix$IGH_Reads[which(clin=="STA")])
mean(summaryMatrix$IGH_Reads[which(clin=="AMR")])
mean(summaryMatrix$IGH_Reads[which(clin=="CMR")])

mean(summaryMatrix$IGK_Reads[which(clin=="STA")])
mean(summaryMatrix$IGK_Reads[which(clin=="AMR")])
mean(summaryMatrix$IGK_Reads[which(clin=="CMR")])

mean(summaryMatrix$IGL_Reads[which(clin=="STA")])
mean(summaryMatrix$IGL_Reads[which(clin=="AMR")])
mean(summaryMatrix$IGL_Reads[which(clin=="CMR")])


COLOR=brewer.pal(3,"Set2")
summaryMatrix<-cbind(data.frame(summaryMatrix),clin)
# Percent TRA & TRB of total T-reads
tiff("Results/MIXCR/Boxplot_IG_expression.tiff",res=300,w=2500,h=2000)
g1<-ggplot(summaryMatrix, aes(clin, IGH_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("IGH expression") + xlab("")
g2<-ggplot(summaryMatrix, aes(clin, IGK_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("IGK expression") + xlab("")
g3<-ggplot(summaryMatrix, aes(clin, IGK_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("IGK expression") + xlab("")
grid.arrange(g1,g2,g3,ncol=2)
dev.off()

####Summary of the TC's
mean(summaryMatrix$TRA_Reads[which(clin=="STA")])
mean(summaryMatrix$TRA_Reads[which(clin=="AMR")])
mean(summaryMatrix$TRA_Reads[which(clin=="CMR")])

mean(summaryMatrix$TRB_Reads[which(clin=="STA")])
mean(summaryMatrix$TRB_Reads[which(clin=="AMR")])
mean(summaryMatrix$TRB_Reads[which(clin=="CMR")])

mean(summaryMatrix$TRG_Reads[which(clin=="STA")])
mean(summaryMatrix$TRG_Reads[which(clin=="AMR")])
mean(summaryMatrix$TRG_Reads[which(clin=="CMR")])

mean(summaryMatrix$TRD_Reads[which(clin=="STA")])
mean(summaryMatrix$TRD_Reads[which(clin=="AMR")])
mean(summaryMatrix$TRD_Reads[which(clin=="CMR")])

tiff("Results/MIXCR/Boxplot_T_expression.tiff",res=300,w=2000,h=2000)
boxplot(summaryMatrix$T_expression~clin,col=COLOR,main="T expression")
dev.off()
summary(glm(summaryMatrix$T_expression~clin)) #p(STA vs AMR)=0.03

tiff("Results/MIXCR/Boxplot_TRs_expression.tiff",res=300,w=2500,h=2000)
g1<-ggplot(summaryMatrix, aes(clin, TRA_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRA expression") + xlab("")
g2<-ggplot(summaryMatrix, aes(clin, TRB_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRB expression") + xlab("")
g3<-ggplot(summaryMatrix, aes(clin, TRG_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRG expression") + xlab("")
g4<-ggplot(summaryMatrix, aes(clin, TRD_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRD expression") + xlab("")
grid.arrange(g1,g2,g3,g4,ncol=2)
dev.off()

summary(glm(summaryMatrix$TRA_expression~relevel(clin,ref="AMR"))) #(STA vs AMR p=0.003, CMR vs AMR p=0.004)
summary(glm(summaryMatrix$TRB_expression~clin)) #(STA vs AMR p=0.04)
summary(glm(summaryMatrix$TRD_expression~clin)) #(STA vs CMR p=0.002 STA vs AMR p=0.0007)
summary(glm(summaryMatrix$TRG_expression~clin)) #(STA vs CMR p=0.1 STA vs AMR p=0.04)

tiff("Results/MIXCR/Boxplot_alphabeta.tiff",res=300,w=2000,h=2000)
ggplot(summaryMatrix, aes(clin, Alpha_Beta_ratio_expression, fill = clin)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab(expression(paste(alpha,beta,"/",gamma,delta,"_ratio"))) + xlab("") + ylim(.85, 1)
dev.off()
summary(glm(summaryMatrix$Alpha_Beta~clin)) #p(STA vs CMR)=0.01 p(STA vs AMR)=0.0001



##### Clinical Data Analysis #######
clin_data<-read.csv("Data/ClinicalData.csv")

id_clin<-match(rownames(summaryMatrix),clin_data$Individual_id)
clin_data<-clin_data[id_clin,]
clin_data$clin<-clin

#donor age
mean(clin_data$donor.age[which(clin_data$clin=="STA")])
mean(clin_data$donor.age[which(clin_data$clin=="CMR")])
mean(clin_data$donor.age[which(clin_data$clin=="AMR")])
sd(clin_data$donor.age[which(clin_data$clin=="STA")])
sd(clin_data$donor.age[which(clin_data$clin=="CMR")])
sd(clin_data$donor.age[which(clin_data$clin=="AMR")])
summary(lm(clin_data$donor.age~clin_data$clin))

#recipient age
mean(clin_data$Rec.age[which(clin_data$clin=="STA")])
mean(clin_data$Rec.age[which(clin_data$clin=="CMR")])
mean(clin_data$Rec.age[which(clin_data$clin=="AMR")])
sd(clin_data$Rec.age[which(clin_data$clin=="STA")])
sd(clin_data$Rec.age[which(clin_data$clin=="CMR")])
sd(clin_data$Rec.age[which(clin_data$clin=="AMR")])
summary(lm(clin_data$Rec.age~clin_data$clin))

#recipient gender
table(clin_data$Rec.Gender,clin_data$clin)
chisq.test(table(clin_data$Rec.Gender,clin_data$clin))

#donor gender
table(clin_data$Donor..gender,clin_data$clin)
chisq.test(table(clin_data$Donor..gender,clin_data$clin))

#Number of TX
mean(clin_data$number.TX[which(clin_data$clin=="STA")])
mean(clin_data$number.TX[which(clin_data$clin=="CMR")])
mean(clin_data$number.TX[which(clin_data$clin=="AMR")])
sd(clin_data$number.TX[which(clin_data$clin=="STA")])
sd(clin_data$number.TX[which(clin_data$clin=="CMR")])
sd(clin_data$number.TX[which(clin_data$clin=="AMR")])
summary(lm(clin_data$number.TX~clin_data$clin))

chisq.test(table(clin_data$Induction..Bxb..1..or.rATg..2..none..0..,clin_data$clin))

chisq.test(table(clin_data$CNI.time.rejection..FK.1.CsA.2.none.0.,clin_data$clin))

chisq.test(table(clin_data$DSA_HLA.class..0.neg..1.class.I..2.class.II..3.class.I.and.II.,clin_data$clin))

chisq.test(table(clin_data$Anti.HLA.Ab_HLA.class,clin_data$clin))

#eGFR
mean(clin_data$eGFR..time.Bx.[which(clin_data$clin=="STA")])
mean(clin_data$eGFR..time.Bx.[which(clin_data$clin=="CMR")])
mean(clin_data$eGFR..time.Bx.[which(clin_data$clin=="AMR")])
sd(clin_data$eGFR..time.Bx.[which(clin_data$clin=="STA")])
sd(clin_data$eGFR..time.Bx.[which(clin_data$clin=="CMR")])
sd(clin_data$eGFR..time.Bx.[which(clin_data$clin=="AMR")])
summary(lm(clin_data$eGFR..time.Bx.~clin_data$clin))

#GraftLoss
chisq.test(table(clin_data$GraftLossCateg,clin_data$clin))


summary(glm(summaryMatrix$Alpha_Beta_ratio_expression~clin_data$GraftLossCateg)) #
boxplot(summaryMatrix$Alpha_Beta~clin_data$donor.age,col=COLOR)



#####################
#### GTEX data ######
####################
load("Data/RepertoireResults_GTEx.Rdata")
COLOR=brewer.pal(8,"Set2")

####Ratio
ratio<-c(summaryMatrix$Alpha_Beta_ratio_expression,summaryMatrix_GTEX$AlphaBeta_Percentage)
clin2<-factor(c(as.character(clin),rep("GTEx",length(summaryMatrix_GTEX$AlphaBeta_Percentage))))
clin2<-factor(clin2,levels=c("STA","CMR","AMR","GTEx"))
summary(glm(ratio~relevel(clin2,ref="AMR"))) 
ratio_TX_GTEX<-cbind(data.frame(ratio),clin2)
tiff("Results/MIXCR/Boxplot_alphabeta_GTEx.tiff",res=300,w=2000,h=2000)
ggplot(ratio_TX_GTEX, aes(clin2, ratio, fill = clin2)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab(expression(paste(alpha,beta,"/",gamma,delta,"_ratio"))) + xlab("")
dev.off()

gtex_id<-read.csv("Data/GTEX/SraRunTable_blood_1691.csv")
names_short<-rownames(summaryMatrix_GTEX)[280:nrow(summaryMatrix_GTEX)]
splitpop2<-strsplit(names_short,"}")
names_short<-unlist(lapply(splitpop2, "[", 1))
rownames(summaryMatrix_GTEX)[280:nrow(summaryMatrix_GTEX)]<-names_short
individual_id<-gtex_id[match(rownames(summaryMatrix_GTEX),gtex_id$Run_s),"submitted_subject_id_s"]
gtex_phenotype_data<-read.csv("Data/GTEX/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.csv")

###clasify diseases by renal, others and unknown
COD<-gtex_phenotype_data[match(individual_id,gtex_phenotype_data$SUBJID),"CauseOfDeath"]
clin3<-factor(c(as.character(clin),as.character(COD)))
clin3<-factor(clin3,levels=c("STA","CMR","AMR","Renal","Liver","Respiratory","Cerebrovascular","Cardio","Trauma","Neuro","Unknown/Other"))
summary(glm(ratio~relevel(clin3,ref="AMR"))) #

###Renal vs Other
COD2<-ifelse(COD=="Renal","GTEx-Renal","GTEx-Other")
clin4<-factor(c(as.character(clin),as.character(COD2)))
clin4<-factor(clin4,levels=c("STA","CMR","AMR","GTEx-Renal","GTEx-Other"))
summary(glm(ratio~relevel(clin4,ref="AMR"))) #
ratio_TX_GTEX<-cbind(data.frame(ratio),clin4)
tiff("Results/MIXCR/Boxplot_alphabeta_GTEx_renal_other.tiff",res=300,w=2000,h=2000)
ggplot(ratio_TX_GTEX, aes(clin4, ratio, fill = clin4)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab(expression(paste(alpha,beta,"/",gamma,delta,"_ratio"))) + xlab("")
dev.off()


####TR's expression
tiff("Results/MIXCR/Boxplot_TRs_expression_GTEx.tiff",res=300,w=3500,h=2500)

TRA<-c(summaryMatrix$TRA_expression,summaryMatrix_GTEX$TRA_expression)
TRA_TX_GTEX<-cbind(data.frame(TRA),clin4)
g1<-ggplot(TRA_TX_GTEX, aes(clin2, TRA, fill = clin2)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRA expression") + xlab("")
summary(glm(TRA~relevel(clin2,ref="GTEx")))

TRB<-c(summaryMatrix$TRB_expression,summaryMatrix_GTEX$TRB_expression)
TRB_TX_GTEX<-cbind(data.frame(TRB),clin2)
g2<-ggplot(TRB_TX_GTEX, aes(clin2, TRB, fill = clin2)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRB expression") + xlab("")
summary(glm(TRB~relevel(clin2,ref="GTEx")))

TRG<-c(summaryMatrix$TRG_expression,summaryMatrix_GTEX$TRG_expression)
TRG_TX_GTEX<-cbind(data.frame(TRG),clin2)
g3<-ggplot(TRG_TX_GTEX, aes(clin2, TRG, fill = clin2)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRG expression") + xlab("")
summary(glm(TRG~relevel(clin2,ref="GTEx")))

TRD<-c(summaryMatrix$TRD_expression,summaryMatrix_GTEX$TRD_expression)
TRD_TX_GTEX<-cbind(data.frame(TRD),clin2)
g4<-ggplot(TRD_TX_GTEX, aes(clin2, TRD, fill = clin2)) + 
  geom_boxplot() + scale_fill_manual(values = COLOR) +
  ylab("TRD expression") + xlab("")
summary(glm(TRD~relevel(clin2,ref="GTEx")))

grid.arrange(g1,g2,g3,g4,ncol=2)
dev.off()

###Summary of the IG's
mean(summaryMatrix_GTEX$IGH_Reads)
mean(summaryMatrix_GTEX$IGK_Reads)
mean(summaryMatrix_GTEX$IGL_Reads)


###Summary of the TR's
mean(summaryMatrix_GTEX$TRA_Reads)
mean(summaryMatrix_GTEX$TRB_Reads)
mean(summaryMatrix_GTEX$TRD_Reads)
mean(summaryMatrix_GTEX$TRG_Reads)

mean(summaryMatrix_GTEX)