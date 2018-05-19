rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: QC analysis
###           
### DESCRIP: Prepate RNAseq data for differential expression and for integration with MIXCR
###         
###
### Author: Silvia Pineda
### Date: October, 2017
############################################################################################
library("RColorBrewer")
library("DESeq2")

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/repertoireResults.Rdata")
load("Data/rsemResults.Rdata")
load("Data/norm_data_rlog_filter.Rdata")

annotation<-read.table("Data/gencode.v23.annotation.txt",sep=";") ##60,498
annotation<-annotation[,-7]
colnames(annotation)<-c("Chr","Start","End","Ensemble_id","type_gene","name")

names(clin)<-colnames(norm_data_rlog)
save(count_rsem_genes,FPKM_rsem_genes,norm_data_rlog,annotation,clin,file="Data/RNAseqProcessedNormalized.Rdata")
