rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Extract summary results from MIXCR for validatiom
###           
### DESCRIP: Analysis with the data from the MIXCR output
###         
###
### Author: Silvia Pineda
### Date: June, 2018
############################################################################################
library(entropy)
library(plyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/dataset_MIXCR_val.Rdata")


samples <- c()

Total_Reads <- c()
IGH_Reads <- c()
IGK_Reads <- c()
IGL_Reads <- c()
TRA_Reads <- c()
TRB_Reads <- c()
TRD_Reads <- c()
TRG_Reads <- c()

sampleNames <- unique(alignmentData_val$Sample) # Array of all sample names
samples_full <-c()
for (i in sampleNames)
{
  # Get sample name
  samples <- append(samples, i)
  
  
  # Get sample alignment data
  sampleSet <- subset(alignmentData_val, Sample == i) # Alignment data for sample i
  sampleGenes <- table(sampleSet$Gene) # Gene frequencies in alignments for sample i
  
  if (dim(sampleGenes)>0){
    samples_full<-append(samples_full, i)
    Total_Reads <- append(Total_Reads, nrow(sampleSet))
    IGH_Reads <- append(IGH_Reads, sampleGenes[["IGH"]])
    IGK_Reads <- append(IGK_Reads, sampleGenes[["IGK"]])
    IGL_Reads <- append(IGL_Reads, sampleGenes[["IGL"]])
    TRA_Reads <- append(TRA_Reads, sampleGenes[["TRA"]])
    TRB_Reads <- append(TRB_Reads, sampleGenes[["TRB"]])
    if(length(grep("TRD",names(sampleGenes)))==0){
      TRD_Reads <- append(TRD_Reads,0)
    } else {
      TRD_Reads <- append(TRD_Reads, sampleGenes[["TRD"]])
    }
    
    TRG_Reads <- append(TRG_Reads, sampleGenes[["TRG"]])
  }
}
# Put everything together into one matrix
summaryMatrix_Everything <- data.frame(Total_Reads, IGH_Reads, IGK_Reads, IGL_Reads, TRA_Reads, TRB_Reads, TRD_Reads, TRG_Reads)

# Add additional columns for ratios
T_Reads <- (summaryMatrix_Everything$TRA_Reads + summaryMatrix_Everything$TRB_Reads + summaryMatrix_Everything$TRD_Reads + summaryMatrix_Everything$TRG_Reads)
summaryMatrix_Everything$AlphaBeta_Percentage <- (summaryMatrix_Everything$TRA_Reads + summaryMatrix_Everything$TRB_Reads)/(T_Reads)
summaryMatrix_Everything$KappaLambda_Ratio <- (summaryMatrix_Everything$IGK_Reads / summaryMatrix_Everything$IGL_Reads)

row.names(summaryMatrix_Everything) <- samples_full

###Patient data
patient_data<-read.csv("Data/Validation/SraRunTable.csv")
id<-match(rownames(summaryMatrix_Everything),patient_data$Run)
patient_data<-patient_data[id,]

###Total reads for normalization  
val_totalReads<-read.table("Data/Validation/total_reads_VAL.txt",sep=";")
id<-match(rownames(summaryMatrix_Everything),as.character(val_totalReads$V1))

summaryMatrix<-summaryMatrix_Everything[which(is.na(id)==F),]
val_totalReads_qc<-val_totalReads[na.omit(id),]
MappedReads<-val_totalReads_qc$V2

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


############################################################
## Calculate the gene expression per each indiviudal gene ##
############################################################
Vgene_counts<-data.matrix(table(alignmentData_val$bestVGene,alignmentData_val$Sample))
Vgene_length<-read.csv("Data/VgeneLength.csv")
id_gene<-match(rownames(Vgene_counts),Vgene_length$Gene.name)
#Vgene_length[id_gene,"length"]

Vgene_expression<-matrix(NA,nrow(Vgene_counts),ncol(Vgene_counts))
for(i in 1:dim(Vgene_counts)[1]){
  norm<-MappedReads*Vgene_length[id_gene[i],"length"]
  Vgene_expression[i,]<-1000*(Vgene_counts[i,]/norm)
}
rownames(Vgene_expression)<-rownames(Vgene_counts)
colnames(Vgene_expression)<-colnames(Vgene_counts)

##Save the data
summaryMatrix_VAL<-summaryMatrix
Vgene_expression_VAL<-Vgene_expression
save(summaryMatrix_VAL, Vgene_expression_VAL,file="Data/RepertoireResults_VAL.Rdata")


