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
library("DESeq2")

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
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type) ##60,498

#filtering by low counts 
##keep genes that have counts per milion values above 0.5 in at least two libraries
dds_filter <-  dds[rowSums(fpm(dds) > 0.5) >= 2] 
dds_filter$type <- relevel(dds_filter$type, ref = "STA")
##normalized with rlog
normalized_rlog <- rlog(dds_filter, blind=F) 
norm_data_rlog<-assay(normalized_rlog)
save(norm_data_rlog,clin,file="Data/norm_data_rlog_filter.Rdata")

###To find which ones are coding and which ones are non coding
id_genes<-match(rownames(dds_filter),annotation$Ensemble_id)
table(annotation[id_genes,"type_gene"])

dds <- DESeq(dds_filter) # N.B. Takes ~ 2 min.
res1 <- results(dds, contrast=c("type","AMR", "STA"), alpha = .05) #length(which(res1$padj<0.05)) 5,399
res2 <- results(dds, contrast=c("type","STA","CMR"), alpha = .05) #length(which(res2$padj<0.05)) 0
res3 <- results(dds, contrast=c("type","AMR","CMR"), alpha = .05) #length(which(res3$padj<0.05)) 2,947

##################
### AMR vs STA ###
##################
allgenes <- rownames(res1)

#overexpressed
overexpressed <- allgenes[which(res1$log2FoldChange > 1 & res1$padj < 0.05)]
id_over<-match(overexpressed,annotation$Ensemble_id)
annotation_overexpressed<-annotation$name[id_over] #1353
write.table(annotation_overexpressed, "Results/RNAseq/overexpressed_AMR_STA_genelist.txt", quote = F, row.names = F)

#underexpressed
underexpressed <- allgenes[which(res1$log2FoldChange < -1 & res1$padj < 0.05)]
id_under<-match(underexpressed,annotation$Ensemble_id)
annotation_underexpressed<-annotation$name[id_under] #258
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

tiff(filename = "heatmap_AMR_STA.tiff", width = 3000, height = 2000,  res = 300)
pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors,show_rownames=F)
dev.off()

#################
### AMR vs CMR ##
#################
allgenes <- rownames(res3)

#overexpressed
overexpressed <- allgenes[which(res3$log2FoldChange > 1 & res3$padj < 0.05)]
id_over<-match(overexpressed,annotation$Ensemble_id)
annotation_overexpressed<-annotation$name[id_over] #480
write.table(annotation_overexpressed, "Results/RNAseq/overexpressed_AMR_CMR_genelist.txt", quote = F, row.names = F)

#underexpressed
underexpressed <- allgenes[which(res3$log2FoldChange < -1 & res3$padj < 0.05)]
id_under<-match(underexpressed,annotation$Ensemble_id)
annotation_underexpressed<-annotation$name[id_under] #105
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

tiff(filename = "heatmap_AMR_CMR.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors,show_rownames=F)
dev.off()



## 2. Using Random Forest to classify the three categories (VSURF to select the gene of interest)
set.seed(33)

fit<-VSURF(x = t(norm_data_rlog), y=clin, parallel = TRUE,ncores=4)
save(fit,file="Results/RNAseq/Genes_VSURF.Rdata")

load("Results/RNAseq/Genes_VSURF.Rdata")
load("Data/norm_data_rlog_filter.Rdata")
genes<-data.frame(norm_data_rlog[fit$varselect.interp,]) ##33 genes

###Plot the results
id_gene<-match(annotation$Ensemble_id,rownames(genes))
significantResults<-genes[na.omit(id_gene),]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR [1] ,"AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "Results/RNAseq/heatmap_RF.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors,show_rownames=F)
dev.off()






