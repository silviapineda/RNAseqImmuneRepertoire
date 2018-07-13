rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### Normalization for the MIXCR summary matrix output
###
### Author: Silvia Pineda
### Date: July, 2018 
############################################################################################

#### The data needs to be normalized by the total number of reads 
load("SummaryMatrixReadsFromMIXCR.Rdata")
totalReads<-read.table("total_reads.txt",sep=";")

##Map the samples
id<-match(rownames(summaryMatrix),as.character(totalReads$V1))
summaryMatrix<-summaryMatrix[which(is.na(id)==F),]
totalReads_qc<-totalReads[na.omit(id),]
MappedReads<-totalReads_qc$V2

##Total IG reads
summaryMatrix$IG_Reads<-summaryMatrix$IGH_Reads+summaryMatrix$IGK_Reads+summaryMatrix$IGL_Reads
summaryMatrix$T_Reads<-summaryMatrix$TRA_Reads+summaryMatrix$TRB_Reads+summaryMatrix$TRD_Reads+summaryMatrix$TRG_Reads

###Analysis in Overall Ab expression
summaryMatrix$IG_expression<-summaryMatrix$IG_Reads/MappedReads
summaryMatrix$IGH_expression<-summaryMatrix$IGH_Reads/MappedReads
summaryMatrix$IGK_expression<-summaryMatrix$IGK_Reads/MappedReads
summaryMatrix$IGL_expression<-summaryMatrix$IGL_Reads/MappedReads

summaryMatrix$T_expression<-summaryMatrix$T_Reads/MappedReads
summaryMatrix$TRA_expression<-summaryMatrix$TRA_Reads/MappedReads
summaryMatrix$TRB_expression<-summaryMatrix$TRB_Reads/MappedReads
summaryMatrix$TRD_expression<-summaryMatrix$TRD_Reads/MappedReads
summaryMatrix$TRG_expression<-summaryMatrix$TRG_Reads/MappedReads

###Ratio
summaryMatrix$Alpha_Beta_ratio_expression<-(summaryMatrix$TRA_expression+summaryMatrix$TRB_expression)/summaryMatrix$T_expression
summaryMatrix$KappaLambda_ratio_expression <- (summaryMatrix$IGK_expression / summaryMatrix$IGL_expression)

save(summaryMatrix,file="RepertoireResults.Rdata")
