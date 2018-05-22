rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Extract summary results from MIXCR 
###           
### DESCRIP: Analysis with the data from the MIXCR output
###         
###
### Author: Akshay Ravoor and Silvia Pineda
### Date: summer, 2017 
############################################################################################
library(entropy)
library(plyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/dataset_MIXCR_GTEX.Rdata")


samples <- c()

Total_Reads <- c()
IGH_Reads <- c()
IGK_Reads <- c()
IGL_Reads <- c()
TRA_Reads <- c()
TRB_Reads <- c()
TRD_Reads <- c()
TRG_Reads <- c()

Total_Clones <- c()
IGH_Clones <- c()
IGK_Clones <- c()
IGL_Clones <- c()
TRA_Clones <- c()
TRB_Clones <- c()
TRD_Clones <- c()
TRG_Clones <- c()

IGH_Entropy <- c()
IGK_Entropy <- c()
IGL_Entropy <- c()
TRA_Entropy <- c()
TRB_Entropy <- c()
TRD_Entropy <- c()
TRG_Entropy <- c()

sampleNames <- unique(alignmentData$Sample) # Array of all sample names

for (i in sampleNames)
{
  # Get sample name
  samples <- append(samples, i)
  print(i)
  
  
  # Get sample alignment data
  sampleSet <- subset(alignmentData, Sample == i) # Alignment data for sample i
  sampleGenes <- table(sampleSet$Gene) # Gene frequencies in alignments for sample i
  
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
  
  
  # Get sample clonal data
  # sampleClones <- subset(clonalData, Sample == i) # Clonal data for sample i
  # 
  # Total_Clones <- append(Total_Clones, sum(sampleClones$cloneCount))
  # IGH_Clones <- append(IGH_Clones, sum(sampleClones[which(sampleClones$Gene == "IGH"), 2]))
  # IGK_Clones <- append(IGK_Clones, sum(sampleClones[which(sampleClones$Gene == "IGK"), 2]))
  # IGL_Clones <- append(IGL_Clones, sum(sampleClones[which(sampleClones$Gene == "IGL"), 2]))
  # TRA_Clones <- append(TRA_Clones, sum(sampleClones[which(sampleClones$Gene == "TRA"), 2]))
  # TRB_Clones <- append(TRB_Clones, sum(sampleClones[which(sampleClones$Gene == "TRB"), 2]))
  # TRD_Clones <- append(TRD_Clones, sum(sampleClones[which(sampleClones$Gene == "TRD"), 2]))
  # TRG_Clones <- append(TRG_Clones, sum(sampleClones[which(sampleClones$Gene == "TRG"), 2]))
  # 
  # 
  # # Get sample entropy data
  # IGH_Entropy <- append(IGH_Entropy, entropy(sampleClones[which(sampleClones$Gene == "IGH"), 2]))
  # IGK_Entropy <- append(IGK_Entropy, entropy(sampleClones[which(sampleClones$Gene == "IGK"), 2]))
  # IGL_Entropy <- append(IGL_Entropy, entropy(sampleClones[which(sampleClones$Gene == "IGL"), 2]))
  # TRA_Entropy <- append(TRA_Entropy, entropy(sampleClones[which(sampleClones$Gene == "TRA"), 2]))
  # TRB_Entropy <- append(TRB_Entropy, entropy(sampleClones[which(sampleClones$Gene == "TRB"), 2]))
  # TRD_Entropy <- append(TRD_Entropy, entropy(sampleClones[which(sampleClones$Gene == "TRD"), 2]))
  # TRG_Entropy <- append(TRG_Entropy, entropy(sampleClones[which(sampleClones$Gene == "TRG"), 2]))
}

# Put everything together into one matrix
summaryMatrix_Everything <- data.frame(Total_Reads, IGH_Reads, IGK_Reads, IGL_Reads, TRA_Reads, TRB_Reads, TRD_Reads, TRG_Reads)

# Add additional columns for ratios
T_Reads <- (summaryMatrix_Everything$TRA_Reads + summaryMatrix_Everything$TRB_Reads + summaryMatrix_Everything$TRD_Reads + summaryMatrix_Everything$TRG_Reads)
summaryMatrix_Everything$AlphaBeta_Percentage <- (summaryMatrix_Everything$TRA_Reads + summaryMatrix_Everything$TRB_Reads)/(T_Reads)
summaryMatrix_Everything$KappaLambda_Ratio <- (summaryMatrix_Everything$IGK_Reads / summaryMatrix_Everything$IGL_Reads)

row.names(summaryMatrix_Everything) <- samples
#names_short<-rownames(summaryMatrix_Everything)[280:nrow(summaryMatrix_Everything)]
#splitpop2<-strsplit(names_short,"}")
#names_short<-unlist(lapply(splitpop2, "[", 1))
#rownames(summaryMatrix_Everything)[280:nrow(summaryMatrix_Everything)]<-names_short

#Save the data
save(summaryMatrix_Everything, file = "Data/SummaryMatrixReadsFromMIXCR_GTEX.Rdata")

#### The data needs to be normalized by the unmapped reads 
load("Data/SummaryMatrixReadsFromMIXCR_GTEX.Rdata")
gtex_totalReads<-read.table("Data/GTEX/total_reads_GTEX.txt",sep=";")
ids<-read.csv("Data/GTEX/SraRunTable_blood_1691_Transcriptomic.csv")

##Find the samples we have in getex_gene_counts
id<-match(rownames(summaryMatrix_Everything),as.character(gtex_totalReads$V1))

summaryMatrix<-summaryMatrix_Everything[which(is.na(id)==F),]
gtex_totalReads_qc<-gtex_totalReads[na.omit(id),]
MappedReads<-gtex_totalReads_qc$V2

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
Vgene_counts<-data.matrix(table(alignmentData$bestVGene,alignmentData$Sample))
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

summaryMatrix_GTEX<-summaryMatrix
Vgene_expression_GTEX<-Vgene_expression
save(summaryMatrix_GTEX, Vgene_expression_GTEX,file="Data/RepertoireResults_GTEx.Rdata")


