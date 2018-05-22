rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Prepare data from MIXCR tool from GTEX data
###           
### DESCRIP: Analysis with the data from the MIXCR GTEX data
###         
###
### Author: Silvia Pineda
### Date: May, 2018
############################################################################################
library(dplyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

########################
### Alignment Matrix ###
########################

alignmentFiles <- list.files("Data/GTEX/MIXCR_alignments/") ##There is a mixed with RNAseq, WGS and WXS

###Find those that are RNAseq
SraRunTable<-read.csv("Data/GTEX/SraRunTable_blood_1691.csv") ##File from the run
table(SraRunTable$Assay_Type_s) ##456 RNAseq

##Extract the RSS ID
splitpop <- strsplit(alignmentFiles,"_")
rs_id<-unlist(lapply(splitpop, "[", 1))
splitpop2<-strsplit(rs_id,"}")
rs_id2<-unlist(lapply(splitpop2, "[", 1))
rs_id<-substring(rs_id2,2)

##Finr the RNA-seq ones
data_mixcr<-SraRunTable[match(rs_id,SraRunTable$Run_s),]
data_mixcr_rnaseq<-data_mixcr[which(data_mixcr$Assay_Type_s=="RNA-Seq"),]

##Read only the RNAseq ones
alignmentFiles<-alignmentFiles[match(data_mixcr_rnaseq$Run_s,rs_id)]

alignmentData <- c()
for (i in alignmentFiles)
{
  print(i)
  
  alignment <- read.delim(paste("Data/GTEX/MIXCR_alignments/", i, sep = ""), quote = "")
  alignment$Sample <- substr(i, 1, nchar(i)-15)
  
  alignmentData <- rbind(alignmentData, alignment)
}

alignGeneV <- substr(alignmentData[, "bestVGene"], 1, 3); alignGeneV[alignGeneV == ""] <- NA
alignGeneD <- substr(alignmentData[, "bestDGene"], 1, 3); alignGeneD[alignGeneD == ""] <- NA
alignGeneJ <- substr(alignmentData[, "bestJGene"], 1, 3); alignGeneJ[alignGeneJ == ""] <- NA

alignmentData$Gene <- coalesce(alignGeneV, alignGeneD, alignGeneJ)

alignmentData$Sample<-substring(alignmentData$Sample,2,11)

# Save the data
save(alignmentData, file = "Data/dataset_MIXCR_GTEX.Rdata")
