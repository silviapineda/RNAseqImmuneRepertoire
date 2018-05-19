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
library("pheatmap")
library("plyr")
library("RColorBrewer")

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/RNAseqProcessedNormalized.Rdata")


########################################
## Differential expression analysis ####
#######################################

###Considering coding and non-coding genes

###################
## 1. Dseq2 #######
###################
coldata<-read.csv("Data/colData.txt")
rownames(coldata)<-coldata[,1]
coldata<-coldata[,-1]
coldata<-coldata[order(rownames(coldata)),]
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type)
#filtering by low counts 
##keep genes that have counts per milion values above 0.5 in at least two libraries
dds_filter <-  dds[rowSums(fpm(dds) > 0.5) >= 2] 
dds$type <- relevel(dds$type, ref = "STA")

dds <- DESeq(dds) # N.B. Takes ~ 2 min.
res1 <- results(dds, contrast=c("type","AMR", "STA"), alpha = .05) #length(which(res1$padj<0.05)) 5,376
res2 <- results(dds, contrast=c("type","STA","CMR"), alpha = .05) #length(which(res2$padj<0.05)) 0
res3 <- results(dds, contrast=c("type","AMR","CMR"), alpha = .05) #length(which(res3$padj<0.05)) 2,920


### AMR vs STA
allgenes <- rownames(res1)

#overexpressed
overexpressed <- allgenes[which(res1$log2FoldChange > 1 & res1$padj < 0.05)]
id_over<-match(overexpressed,annotation$Ensemble_id)
annotation_overexpressed<-annotation$name[id_over] #1342
write.table(annotation_overexpressed, "Results/RNAseq/overexpressed_AMR_STA_genelist.txt", quote = F, row.names = F)

#underexpressed
underexpressed <- allgenes[which(res1$log2FoldChange < -1 & res1$padj < 0.05)]
id_under<-match(underexpressed,annotation$Ensemble_id)
annotation_underexpressed<-annotation$name[id_under] #259
write.table(annotation_underexpressed, "Results/RNAseq/underexpressed_AMR_STA_genelist.txt", quote = F, row.names = F)


###Plot the results
significantResults <- norm_data_rlog[match(c(overexpressed,underexpressed),rownames(norm_data_rlog)),which(clin!="CMR")]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin[which(clin!="CMR")])
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR[1], "AMR" = COLOR[3]))

tiff(filename = "heatmap_AMR_STA.tiff", width = 15, height = 7, units = "in", res = 300)
pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors)
dev.off()

### AMR vs CMR
allgenes <- rownames(res3)

#overexpressed
overexpressed <- allgenes[which(res3$log2FoldChange > 1 & res3$padj < 0.05)]
id_over<-match(overexpressed,annotation$Ensemble_id)
annotation_overexpressed<-annotation$name[id_over] #476
write.table(annotation_overexpressed, "Results/RNAseq/overexpressed_AMR_CMR_genelist.txt", quote = F, row.names = F)

#underexpressed
underexpressed <- allgenes[which(res3$log2FoldChange < -1 & res3$padj < 0.05)]
id_under<-match(underexpressed,annotation$Ensemble_id)
annotation_underexpressed<-annotation$name[id_under] #106
write.table(annotation_underexpressed, "Results/RNAseq/underexpressed_AMR_CMR_genelist.txt", quote = F, row.names = F)

###Plot the results
significantResults <- norm_data_rlog[match(c(overexpressed,underexpressed),rownames(norm_data_rlog)),which(clin!="STA")]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin[which(clin!="STA")])
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "heatmap_AMR_CMR.tiff", width = 15, height = 7, units = "in", res = 300)
pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors)
dev.off()



## 1. Using Random Forest to classify the three categories (VSURF to select the gene of interest)
set.seed(112233)

fit<-VSURF(x = t(norm_data_rlog_coding), y=clin, parallel = TRUE,ncores=4)
save(fit,file="Results/RNAseq/CodingGenes_RF.Rdata")

genes<-data.frame(norm_data_rlog_coding[fit$varselect.interp,])

set.seed(1000)
rf_output <- randomForest(clin~.,data=t(genes),proximity=TRUE, keep.forest=T,ntree=1000)



#### Coding genes 
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id_gene<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_gene),]







