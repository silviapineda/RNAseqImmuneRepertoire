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
### DESCRIP: Analysis of the output from the alignment
###         
###
### Author: Silvia Pineda
### Date: October, 2017
############################################################################################
library("RColorBrewer")
library("DESeq2")
library(plyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

capture<-read.table("QC/capture.txt",sep="\t")
##Delete the Mix
capture<-capture[which(capture$V3!="M1Log.final.out" & capture$V3!="M2Log.final.out" &
                       capture$V3!="M3Log.final.out" & capture$V3!="M4Log.final.out"),]

library(reshape2)
capture.df<-unstack(capture, form=V2~V1)
rownames(capture.df)<-unique(capture$V3)
colnames(capture.df)<-c("numReads","percReadsUnmappedOther","uniquelyMappedReads","percReadsUnmappedTooShort","percReadsUnmappedManyMismatches")
rownames(capture.df)<-gsub('.{13}$', '', rownames(capture.df))


clin<-factor(substr(rownames(capture.df),1,1))
clin<-revalue(clin, c("C"="CMR", "H"="AMR" , "N" = "STA"))
clin<-relevel(clin,ref="STA")
COLOR=brewer.pal(3,"Set2")
cols = COLOR[clin]
tiff("NumberReads.tiff",res=300,w=4000,h=2000)
barplot(as.numeric(capture.df$numRead[order(as.numeric(capture.df$numReads))]),
        col = cols[order(as.numeric(capture.df$numReads))],
        names.arg = rownames(capture.df)[order(as.numeric(capture.df$numReads))],
        cex.names=0.5,las=2)
legend(0, 80000000, legend=levels(clin),col=COLOR,pch=15, cex=1)
dev.off()

tiff("UniqueMappedReads.tiff",res=300,w=4000,h=2000)
barplot(as.numeric(capture.df$uniquelyMappedReads[order(as.numeric(capture.df$uniquelyMappedReads))]),
        col = cols[order(as.numeric(capture.df$uniquelyMappedReads))],
        names.arg = rownames(capture.df)[order(as.numeric(capture.df$uniquelyMappedReads))],
        cex.names=0.5,las=2)
legend(0, 70000000, legend=levels(clin),col=COLOR,pch=15, cex=1)
dev.off()

capture.df$percReadsUnmappedTooShort <- as.numeric(sub("%", "", capture.df$percReadsUnmappedTooShort))/100
capture.df$percReadsUnmappedOther <- as.numeric(sub("%", "", capture.df$percReadsUnmappedOther))/100

capture.df$numUnmappedReads<-as.numeric(capture.df$numReads)*(capture.df$percReadsUnmappedTooShort+capture.df$percReadsUnmappedOther)

#save(capture.df,file="~/ImmuneRep_RNAseq/MIXCR/summaryCapture.Rdata")

tiff("QC/UnmappeddReads.tiff",res=300,w=4000,h=2000)
barplot(as.numeric(capture.df$numUnmappedReads[order(as.numeric(capture.df$numUnmappedReads))]),
        col = cols[order(as.numeric(capture.df$numUnmappedReads))],
        names.arg = rownames(capture.df)[order(as.numeric(capture.df$numUnmappedReads))],
        cex.names=0.5,las=2)
legend(0, 3000000, legend=levels(clin),col=COLOR,pch=15, cex=1)
dev.off()

####Quality control considering the MIXCR output
load("MIXCR/summaryMatrix.Rdata")
id<-match(rownames(summaryMatrix),rownames(capture.df))

tiff("Plot_Reads_UnmappeddReads.tiff",res=300,w=2000,h=2000)
plot(summaryMatrix$Total_Reads,capture.df$numUnmappedReads[id],col=cols,pch=19)
cor(summaryMatrix$Total_Reads,capture.df$numUnmappedReads[id]) #0.8
legend(1500, 3300000, legend=levels(clin),col=COLOR,pch=19, cex=1)
dev.off()

tiff("Plot_Ratio_UnmappeddReads.tiff",res=300,w=2000,h=2000)
plot(summaryMatrix$Alpha_Beta,capture.df$numUnmappedReads[id],col=cols,pch=19)
cor(summaryMatrix$Alpha_Beta,capture.df$numUnmappedReads[id]) #0.2
legend(0.87, 3300000, legend=levels(clin),col=COLOR,pch=19, cex=1)
dev.off()

tiff("Plot_Reads_NumdReads.tiff",res=300,w=2000,h=2000)
plot(summaryMatrix$Total_Reads,capture.df$numReads[id],col=cols,pch=19)
cor(summaryMatrix$Total_Reads,as.numeric(capture.df$numReads[id])) #0.8
legend(1500, 85000000, legend=levels(clin),col=COLOR,pch=19, cex=1)
dev.off()

tiff("Plot_Ratio_NumdReads.tiff",res=300,w=2000,h=2000)
plot(summaryMatrix$Alpha_Beta,capture.df$numReads[id],col=cols,pch=19)
cor(summaryMatrix$Alpha_Beta,as.numeric(capture.df$numReads[id])) #0.8
legend(0.87, 85000000, legend=levels(clin),col=COLOR,pch=19, cex=1)
dev.off()


####### QC3 in fastq files
qc3_report<-read.table("QC/QC3/fastqSummary2.txt",header=T)
summary(qc3_report$TotalReads)
summary(qc3_report$GC)
summary(qc3_report$BQ)

clin<-factor(substr(qc3_report$Sample_ID,1,1))
clin<-revalue(clin, c("C"="CMR", "H"="AMR" , "N" = "STA"))
clin<-relevel(clin,ref="STA")
COLOR=brewer.pal(3,"Set2")
cols = COLOR[clin]
tiff("NumberReads_fastq.tiff",res=300,w=4000,h=2000)
barplot(qc3_report$TotalReads[order(qc3_report$TotalReads)],
        col = cols[order(qc3_report$TotalReads)],
        names.arg = qc3_report$Sample_ID[order(qc3_report$TotalReads)],
        cex.names=0.5,las=2)
legend(0, 80000000, legend=levels(clin),col=COLOR,pch=15, cex=1)
dev.off()


####### RSEM ###########
files <- list.files("RSEM/RSEMnoMixed/")

count_rsem_genes <- c()
TPM_rsem_genes <- c()
FPKM_rsem_genes <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("RSEM/RSEMnoMixed/", i, sep = ""))
  count_rsem_genes <- cbind(count_rsem_genes, t[,5])
  TPM_rsem_genes <- cbind(TPM_rsem_genes, t[,6])
  FPKM_rsem_genes <- cbind(FPKM_rsem_genes, t[,7])
}

rownames(count_rsem_genes)<-t[,1]
colnames(count_rsem_genes)<-gsub('.{18}$', '', files)

rownames(TPM_rsem_genes)<-t[,1]
colnames(TPM_rsem_genes)<-gsub('.{18}$', '', files)

rownames(FPKM_rsem_genes)<-t[,1]
colnames(FPKM_rsem_genes)<-gsub('.{18}$', '', files)

###Delete the miixed samples
#count_rsem_genes<-count_rsem_genes[,-c(which(colnames(count_rsem_genes)=="M2" 
#                                         | colnames(count_rsem_genes)=="M3" | colnames(count_rsem_genes)=="M4" ))]
#TPM_rsem_genes<-TPM_rsem_genes[,-c(which(colnames(TPM_rsem_genes)=="M2" 
#                                             | colnames(TPM_rsem_genes)=="M3" | colnames(TPM_rsem_genes)=="M4" ))]
#FPKM_rsem_genes<-FPKM_rsem_genes[,-c(which(colnames(FPKM_rsem_genes)=="M2" 
#                                             | colnames(FPKM_rsem_genes)=="M3" | colnames(FPKM_rsem_genes)=="M4" ))]


##Save them ordered to make them comparable
count_rsem_genes<-count_rsem_genes[,order(colnames(count_rsem_genes))]
TPM_rsem_genes<-TPM_rsem_genes[,order(colnames(TPM_rsem_genes))]
FPKM_rsem_genes<-FPKM_rsem_genes[,order(colnames(FPKM_rsem_genes))]
save(count_rsem_genes,TPM_rsem_genes,FPKM_rsem_genes, file="rsemResults.Rdata")

#####


#########
## Some QC in the counts and normalization 
#########
load("Data/rsemResults.Rdata")

#To obtain the PCA
pca<-prcomp(count_rsem_genes)
clin<-factor(substr(colnames(count_rsem_genes),1,1))
clin<-revalue(clin, c("C"="CMR", "H"="AMR" , "N" = "STA"))
clin<-relevel(clin,ref="STA")
SPP <- clin
COLOR=brewer.pal(3,"Set2") 
tiff("PCAplot_counts.tiff",res=300,w=2000,h=2000)
plot(pca$rotation[,1],pca$rotation[,2],col=COLOR[SPP],pch=20,cex=1.5,xlab="PC1",ylab="PC2")
dev.off()

pca<-prcomp(FPKM_rsem_genes_filter)
clin<-factor(substr(colnames(FPKM_rsem_genes),1,1))
clin<-revalue(clin, c("C"="CMR", "H"="AMR" , "N" = "STA"))
clin<-relevel(clin,ref="STA")
SPP <- clin
COLOR=brewer.pal(3,"Set2") 
tiff("PCAplot_FPKM.tiff",res=300,w=2000,h=2000)
plot(pca$rotation[,1],pca$rotation[,2],col=COLOR[SPP],pch=20,cex=1.5,xlab="PC1",ylab="PC2")
dev.off()

##Batch effects
batch_effects<-read.table("/Users/Pinedasans/ImmuneRep_RNAseq/QC/QC3/BatchEffect.txt",header=T)
batch_effects<-batch_effects[order(batch_effects$Sample_ID),]
SPP<-batch_effects$Instrument
COLOR=brewer.pal(3,"Set1") 
tiff("PCAplot_FPKM_batchEffect.tiff",res=300,w=2000,h=2000)
plot(pca$rotation[,1],pca$rotation[,2],col=COLOR[SPP],pch=20,xlab="PC1",ylab="PC2",cex=1.5)
dev.off()


###Boxplots for nomalizations
FPKM_rsem_genes_filter <- FPKM_rsem_genes[ abs(rowSums(FPKM_rsem_genes)) > 0, ] 
TPM_rsem_genes_filter <- TPM_rsem_genes[ abs(rowSums(TPM_rsem_genes)) > 0, ] 
tiff("Boxplot_FPKM.tiff",res=300,w=4000,h=2000)
boxplot(log2(FPKM_rsem_genes_filter+0.001),col=COLOR[SPP],cex.axis=0.6)
dev.off()
tiff("Boxplot_TPM.tiff",res=300,w=4000,h=2000)
boxplot(log2(TPM_rsem_genes_filter+0.001),col=COLOR[SPP],cex.axis=0.6)
dev.off()

##Aplying DSEq to the counts
coldata<-read.csv("RNAseq/colData.txt")
rownames(coldata)<-coldata[,1]
coldata<-coldata[,-1]
coldata<-coldata[order(rownames(coldata)),]
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type)
#filtering by low counts 
##keep genes that have counts per milion values above 0.5 in at least two libraries
dds_filter <-  dds[rowSums(fpm(dds) > 0.5) >= 2] 
##normalized with rlog
normalized_rlog <- rlog(dds_filter, blind=F) 
norm_data_rlog<-assay(normalized_rlog)
save(norm_data_rlog,file="Data/norm_data_filter_rlog_filter.Rdata")

load("Data/norm_data_rlog.Rdata")
tiff("Boxplot_rlog.tiff",res=300,w=4000,h=2000)
boxplot(norm_data_rlog,col=COLOR[SPP],cex.axis=0.6)
dev.off()

pca<-prcomp(norm_data_rlog)
tiff("PCAplot_rlog.tiff",res=300,w=2000,h=2000)
plot(pca$rotation[,1],pca$rotation[,2],col=COLOR[SPP],pch=20,cex=1.5,xlab="PC1",ylab="PC2")
dev.off()


