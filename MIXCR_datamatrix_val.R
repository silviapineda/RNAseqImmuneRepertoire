rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Prepare data from MIXCR tool for the validation datset
###           
### DESCRIP: Analysis with the data from the MIXCR output
###         
###
### Author: Silvia Pineda
### Date: June, 2018
############################################################################################
library(dplyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

########################
### Alignment Matrix ###
########################

alignmentFiles <- list.files("Data/Validation/MIXCR_alignments/")
alignmentData <- c()

for (i in alignmentFiles)
{
  print(i)
  
  alignment <- read.delim(paste("Data/Validation/MIXCR_alignments/", i, sep = ""), quote = "")
  alignment$Sample <- substr(i, 1, nchar(i)-15)
  
  alignmentData <- rbind(alignmentData, alignment)
}

alignGeneV <- substr(alignmentData[, "bestVGene"], 1, 3); alignGeneV[alignGeneV == ""] <- NA
alignGeneD <- substr(alignmentData[, "bestDGene"], 1, 3); alignGeneD[alignGeneD == ""] <- NA
alignGeneJ <- substr(alignmentData[, "bestJGene"], 1, 3); alignGeneJ[alignGeneJ == ""] <- NA

alignmentData$Gene <- coalesce(alignGeneV, alignGeneD, alignGeneJ)

# Save the data
alignmentData_val<-alignmentData
save(alignmentData_val, file = "Data/dataset_MIXCR_val.Rdata")

