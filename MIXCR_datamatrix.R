rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Prepare data from MIXCR tool
###           
### DESCRIP: Analysis with the data from the MIXCR output
###         
###
### Author: Akshay Ravoor
### Date: summer, 2017
############################################################################################
library(dplyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

########################
### Alignment Matrix ###
########################

alignmentFiles <- list.files("Data/MIXCR_alignments/")
alignmentData <- c()

for (i in alignmentFiles)
{
  #print(i)
  
  alignment <- read.delim(paste("Data/MIXCR_alignments/", i, sep = ""), quote = "")
  alignment$Sample <- substr(i, 1, nchar(i)-15)
  
  alignmentData <- rbind(alignmentData, alignment)
}

alignGeneV <- substr(alignmentData[, "bestVGene"], 1, 3); alignGeneV[alignGeneV == ""] <- NA
alignGeneD <- substr(alignmentData[, "bestDGene"], 1, 3); alignGeneD[alignGeneD == ""] <- NA
alignGeneJ <- substr(alignmentData[, "bestJGene"], 1, 3); alignGeneJ[alignGeneJ == ""] <- NA

alignmentData$Gene <- coalesce(alignGeneV, alignGeneD, alignGeneJ)



#####################
### Clonal Matrix ###
#####################

cloneFiles <- list.files("Data/MIXCR_clones/")
clonalData <- c()

for (i in cloneFiles)
{
  #print(i)
  
  clone <- read.delim(paste("Data/MIXCR_clones/", i, sep = ""), quote = "")
  clone$Sample <- substr(i, 1, nchar(i)-11)
  
  clonalData <- rbind(clonalData, clone)
}

cloneGeneV <- substr(clonalData[, "bestVGene"], 1, 3); cloneGeneV[cloneGeneV == ""] <- NA
cloneGeneD <- substr(clonalData[, "bestDGene"], 1, 3); cloneGeneD[cloneGeneD == ""] <- NA
cloneGeneJ <- substr(clonalData[, "bestJGene"], 1, 3); cloneGeneJ[cloneGeneJ == ""] <- NA

clonalData$Gene <- coalesce(cloneGeneV, cloneGeneD, cloneGeneJ)


# Save the data
save(alignmentData, clonalData, file = "Data/dataset_MIXCR.Rdata")
