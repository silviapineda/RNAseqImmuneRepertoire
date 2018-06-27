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
library("VSURF")
library("randomForest")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("DESeq2")
library("glmnet")

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

save(norm_data_rlog,clin,annotation,file="Data/norm_data_rlog_filter.Rdata")

load("Data/norm_data_rlog_filter.Rdata")
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
#1. DSEQ2
allgenes <- rownames(res1)

#results
genes_sign <- allgenes[which(res1$padj < 0.05)]
id_sign<-match(genes_sign,annotation$Ensemble_id)
annotation_sign<-annotation[id_sign,] 
annotation_sign$log2FC<-res1$log2FoldChange[which(res1$padj < 0.05)]
results_STA_AMR_DEseq<-annotation_sign
write.csv(results_STA_AMR_DEseq, "Results/RNAseq/results_STA_AMR_DEseq.csv", row.names = F)

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_Dseq2.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "STA", results_STA_AMR_DEseq, color=COLOR[c(1,3)])
dev.off()

#2. ENET with binomial distribution
results_STA_AMR_ENET<-ENET_binomial(clin, "AMR", "STA", norm_data_rlog) ##60 genes
write.csv(results_STA_AMR_ENET,"Results/RNAseq/results_STA_AMR_ENET.csv",row.names = F)

COLOR = brewer.pal(4,"Set2")
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_ENET.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "STA", results_STA_AMR, color=COLOR[c(1,3)])
dev.off()

id_overlap<-match(results_STA_AMR_ENET$name,results_STA_AMR_DEseq$name)


##################
### CMR vs STA ###
##################
#1. DSEQ2
allgenes <- rownames(res2)

#results
genes_sign <- allgenes[which(res2$padj < 0.05)] ##No genes significant

#2. ENET with binomial distribution
results_STA_CMR_ENET<-ENET_binomial(clin, "CMR", "STA", norm_data_rlog) ## 1 gene
write.csv(results_STA_CMR_ENET,"Results/RNAseq/results_STA_CMR_ENET.csv",row.names = F)

id<-match(results_STA_CMR_ENET$Ensemble_id,rownames(norm_data_rlog))
tiff(filename = "Results/RNAseq/heatmap_CMR_STA_ENET.tiff", width = 2000, height = 2000,  res = 300)
boxplot(norm_data_rlog[id,which(clin!="AMR")]~factor(clin[which(clin!="AMR")]),col=COLOR[c(1,2)])
dev.off()

##################
### AMR vs CMR ###
##################
#1. DSEQ2
allgenes <- rownames(res3)

#results
genes_sign <- allgenes[which(res3$padj < 0.05)]
id_sign<-match(genes_sign,annotation$Ensemble_id)
annotation_sign<-annotation[id_sign,] 
annotation_sign$log2FC<-res3$log2FoldChange[which(res3$padj < 0.05)]
results_CMR_AMR_DEseq<-annotation_sign
write.csv(results_CMR_AMR_DEseq, "Results/RNAseq/results_CMR_AMR_DEseq.csv", row.names = F)

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_CMR_AMR_Dseq2.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "CMR", results_STA_AMR_DEseq, color=COLOR[c(2,3)])
dev.off()

#2. ENET with binomial distribution
results_CMR_AMR_ENET<-ENET_binomial(clin, "AMR", "CMR", norm_data_rlog) ##21 genes
write.csv(results_CMR_AMR_ENET,"Results/RNAseq/results_CMR_AMR_ENET.csv",row.names = F)

COLOR = brewer.pal(4,"Set2")
tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_ENET.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "CMR", results_CMR_AMR_ENET, color=COLOR[c(2,3)])
dev.off()

id_overlap<-match(results_CMR_AMR_ENET$name,results_CMR_AMR_DEseq$name)


#########################
## AMR vs CMR vs STA ###
########################
## 1. Using Random Forest to classify the three categories (VSURF to select the gene of interest)
set.seed(54)
load("Data/norm_data_rlog_filter.Rdata")

fit<-VSURF(x = t(norm_data_rlog), y=clin, parallel = TRUE,ncores=4)
save(fit,file="Results/RNAseq/Genes_VSURF.Rdata")

load("Results/RNAseq/Genes_VSURF.Rdata")

genes<-data.frame(norm_data_rlog[fit$varselect.interp,]) ##10 genes

id<-match(rownames(genes),annotation$Ensemble_id)
rownames(genes)<-annotation$name[id]

###Plot the results
id_gene<-match(annotation$name,rownames(genes))
significantResults<-genes[na.omit(id_gene),]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR [1] ,"AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "Results/RNAseq/heatmap_RF_NoncodingGenes.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()


## 2. Using multinomial ENET
###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(33)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog),clin,family="multinomial",type.multinomial = "grouped"
                                                        ,standardize=TRUE,alpha=a,nfolds=5))})
xx<-rep(NA,length(alphalist))
yy<-rep(NA,length(alphalist))
for (j in 1:length(alphalist)) {
  #print(j)
  if(class(elasticnet[[j]]) != "try-error"){
    xx[j]<-elasticnet[[j]]$lambda.min
    id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
    yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
  }
}
id.min<-which(yy==min(yy,na.rm=TRUE))
lambda<-xx[id.min]
alpha<-alphalist[id.min]

enet<-glmnet(t(norm_data_rlog),clin,family="multinomial",type.multinomial = "grouped",standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta[[2]])[which(enet$beta[[2]]!=0)]
coef1<-enet$beta[[1]][which(enet$beta[[1]]!=0)]
coef2<-enet$beta[[2]][which(enet$beta[[2]]!=0)]
results<-annotation[match(genes,annotation$Ensemble_id),]
write.csv(cbind(results,coef1,coef2),"Results/RNAseq//genes.enet.multinomial.csv",row.names = F)

###Plot the results
id_gene<-match(genes,rownames(norm_data_rlog))
significantResults<-norm_data_rlog[na.omit(id_gene),]

id<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR [1] ,"AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "Results/RNAseq/heatmap_ENET_multinomial_ALL.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()

