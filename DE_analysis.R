rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Diff Expression analysis uisng the RNASeq data
###           
### DESCRIP: Analysis with the data from the QC process
###         
###
### Author: Silvia Pineda
### Date: October, 2017
############################################################################################
library(VSURF)
library(randomForest)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/RNAseqProcessedNormalized.Rdata")

###Read the clinical data and prepare it to match with  the rest of the data
clinical_data<-read.csv("Data/ClinicalData.csv")
colnames(norm_data_rlog)
id_merge<-match(names(clin),clinical_data$Individual_id)
clinical_data<-clinical_data[na.omit(id_merge),]
clinical_data$clin<-clin
clinical_data$Individual_id<-factor(clinical_data$Individual_id)

#######################
## Demographic table ##
#######################



########################################
## Differential expression analysis ####
#######################################

#### Coding genes 
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id_gene<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_gene),]

## 1. Using Random Forest to classify the three categories (VSURF to select the gene of interest)
set.seed(112233)

fit<-VSURF(x = t(norm_data_rlog_coding), y=clin, parallel = TRUE,ncores=4)
save(fit,file="Results/RNAseq/CodingGenes_RF.Rdata")

genes<-data.frame(norm_data_rlog_coding[fit$varselect.interp,])

set.seed(1000)
rf_output <- randomForest(clin~.,data=t(genes),proximity=TRUE, keep.forest=T,ntree=1000)










