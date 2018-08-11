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
library("GISTools")
library(factoextra)
require("cluster")

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
coldata_rej<-read.csv("Data/colDataRej.txt")

rownames(coldata)<-coldata[,1]
coldata<-coldata[,-1]
coldata<-coldata[order(rownames(coldata)),]
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type) ##60,498

rownames(coldata_rej)<-coldata_rej[,1]
coldata<-coldata_rej[,-1]
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
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id_genes<-match(annotation_coding$Ensemble_id,rownames(dds_filter))
dds_filter_coding<-dds_filter[na.omit(id_genes),] #only coding

dds <- DESeq(dds_filter) # N.B. Takes ~ 2 min.
res1 <- results(dds, contrast=c("type","AMR", "STA"), alpha = .05) #length(which(res1$padj<0.05)) 5,399
res2 <- results(dds, contrast=c("type","CMR","STA"), alpha = .05) #length(which(res2$padj<0.05)) 0
res3 <- results(dds, contrast=c("type","AMR","CMR"), alpha = .05) #length(which(res3$padj<0.05)) 2,947

dds<-DESeq(dds_filter_coding)
res4 <- results(dds, contrast=c("type","REJ","STA"), alpha = .05) #length(which(res3$padj<0.05)) 2,947

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

COLOR = brewer.pal(4,"Pastel1")

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_Dseq2.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "STA", results_STA_AMR_DEseq, color=COLOR[c(3,1)])
dev.off()

#2. ENET with binomial distribution
results_STA_AMR_ENET<-ENET_binomial(clin, "AMR", "STA", norm_data_rlog) ##60 genes
write.csv(results_STA_AMR_ENET,"Results/RNAseq/results_STA_AMR_ENET.csv",row.names = F)

##load the results
results_STA_AMR_ENET<-read.csv("Results/RNAseq/results_STA_AMR_ENET.csv")
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_ENET.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "STA", results_STA_AMR_ENET, color=COLOR[c(3,1)])
dev.off()

tiff(filename = "Results/RNAseq/heatmap_AMR_STA_ENET.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "STA", results_STA_AMR_ENET[which(results_STA_AMR_ENET$type_gene=="protein_coding"),], color=COLOR[c(3,1)])
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
tiff(filename = "Results/RNAseq/heatmap_CMR_STA_ENET.tiff", width = 1500, height = 2000,  res = 300)
boxplot(norm_data_rlog[id,which(clin!="AMR")]~factor(clin[which(clin!="AMR")]),col=COLOR[c(3,2)],ylab="TCAF1 log2 gene expression")
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
plot_heatmap(clin,norm_data_rlog, "AMR", "CMR", results_CMR_AMR_DEseq, color=COLOR[c(2,1)])
dev.off()

#2. ENET with binomial distribution
results_CMR_AMR_ENET<-ENET_binomial(clin, "AMR", "CMR", norm_data_rlog) ##21 genes
write.csv(results_CMR_AMR_ENET,"Results/RNAseq/results_CMR_AMR_ENET.csv",row.names = F)

tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_ENET.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "CMR", results_CMR_AMR_ENET, color=COLOR[c(2,1)])
dev.off()

id_overlap<-match(results_CMR_AMR_ENET$name,results_CMR_AMR_DEseq$name)

##################
### Rej vs. STA ###
##################
#results
allgenes <- rownames(res4)
genes_sign <- allgenes[which(res4$padj < 0.05)]
id_sign<-match(genes_sign,annotation$Ensemble_id)
annotation_sign<-annotation[id_sign,] 
annotation_sign$log2FC<-res1$log2FoldChange[which(res4$padj < 0.05)]
results_STA_REJ_DEseq<-annotation_sign #1105
write.csv(results_STA_REJ_DEseq, "Results/RNAseq/results_STA_REJ_DEseq.csv", row.names = F)

COLOR = brewer.pal(4,"Pastel1")

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_REJ_STA_Dseq2.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin2,norm_data_rlog_coding, "REJ", "STA", results_STA_REJ_DEseq, color=COLOR[c(3,1)])
dev.off()

#2. ENET with binomial distribution
clin2<-ifelse(clin=="AMR" | clin=="CMR", "REJ","STA")
#Only in coding genes to compare with microarray
id_coding<-match(annotation$Ensemble_id[which(annotation$type_gene=="protein_coding")],rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_coding),]

results_REJ_STA_ENET<-ENET_binomial(clin2, "REJ", "STA", norm_data_rlog_coding) ##327 genes
write.csv(results_REJ_STA_ENET,"Results/RNAseq/results_REJ_STA_ENE_codingT.csv",row.names = F)

tiff(filename = "Results/RNAseq/heatmap_REJ_STA_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
out<-plot_heatmap(clin2,norm_data_rlog_coding, "REJ", "STA", results_REJ_STA_ENET, color=COLOR[c(2,3)])
dev.off()

###cluster and similarity
results_REJ_STA_ENET<-read.csv("Results/RNAseq/results_REJ_STA_ENE_codingT.csv")
significantResults <- norm_data_rlog_coding[match(results_REJ_STA_ENET$Ensemble_id,rownames(norm_data_rlog_coding)),]
id_gene<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id_gene]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)

hc.cut <- hcut(t(xtst), k = 2, hc_method = "complete",habillage=clin2)
library("clusteval")
cluster_similarity(hc.cut$cluster,clin2) ##1

cluster<-sort(cutree(out$tree_row, k=2))
id_clust<-match(results_REJ_STA_ENET$name,names(cluster))    
results_REJ_STA_ENET$cluster<-cluster[id_clust]
write.csv(results_REJ_STA_ENET,"Results/RNAseq/results_REJ_STA_ENE_codingT.csv",row.names = F)

#########################
## AMR vs CMR vs STA ###
########################
## 1. Using Random Forest to classify the three categories (VSURF to select the gene of interest)

#fit<-VSURF(x = t(norm_data_rlog), y=clin, parallel = TRUE,ncores=4)
#save(fit,file="Results/RNAseq/Genes_VSURF.Rdata")

load("Results/RNAseq/Genes_VSURF.Rdata")

genes<-data.frame(norm_data_rlog[fit$varselect.interp,]) ##33 genes

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
ann_colors = list(Type = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]))

tiff(filename = "Results/RNAseq/heatmap_3categ_RF.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()


## 2. Using multinomial ENET
###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
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

write.csv(cbind(results,coef1,coef2),"Results/RNAseq/results_AMR_CMR_STA_ENET.csv",row.names = F)

###Plot the results
results<-read.csv("Results/RNAseq/results_AMR_CMR_STA_ENET.csv")

id_gene<-match(results$Ensemble_id,rownames(norm_data_rlog))
significantResults<-norm_data_rlog[na.omit(id_gene),]

id<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Pastel1")
ann_colors = list(Type = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]))

tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_STA_ENET.tiff", width = 3000, height = 3500, res = 300)
out<-pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()

cluster<-sort(cutree(out$tree_row, k=3))
id_clust<-match(results$name,names(cluster))    
results$cluster<-cluster[id_clust]
plot(out$tree_row)
abline(h=9, col="red", lty=2, lwd=2)


## Obtain the clusterization
#COLOR = add.alpha(brewer.pal(4,"Set1"),0.3)
COLOR = brewer.pal(4,"Pastel1")
hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete",habillage=clin)
library("clusteval")
cluster_similarity(hc.cut$cluster,clin) ##1
##clsuter plot
tiff(filename = "Results/RNAseq/cluster_ENET_multinomial.tiff", width = 2000, height = 2000, res = 300)
fviz_cluster(hc.cut, t(significantResults),ellipse.type = "norm",habillage=clin,palette=COLOR[c(2,1,3)]) + theme_minimal()
dev.off()
# # Visualize dendrogram
fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE,habillage=clin,palette=COLOR[c(2,1,3)])
# # Visualize silhouhette information
tiff(filename = "Results/RNAseq/silhouette_ENET_multinomial.tiff", width = 2000, height = 2500, res = 300)
fviz_silhouette(hc.cut,palette=COLOR[c(2,1,3)])
dev.off()


###Obtain the FC
Log2FC_AMR_STA<-res1$log2FoldChange
names(Log2FC_AMR_STA)<-rownames(res1)
Log2FC_CMR_STA<-res2$log2FoldChange
names(Log2FC_CMR_STA)<-rownames(res2)
Log2FC_AMR_CMR<-res3$log2FoldChange
names(Log2FC_AMR_CMR)<-rownames(res3)

id_FC<-match(results$Ensemble_id,names(Log2FC_AMR_STA))
results$Log2FC_AMR_STA<-Log2FC_AMR_STA[id_FC]
id_FC<-match(results$Ensemble_id,names(Log2FC_CMR_STA))
results$Log2FC_CMR_STA<-Log2FC_CMR_STA[id_FC]
id_FC<-match(results$Ensemble_id,names(Log2FC_AMR_CMR))
results$Log2FC_AMR_CMR<-Log2FC_AMR_CMR[id_FC]

write.csv(results,"Results/RNAseq/results_AMR_CMR_STA_ENET.csv",row.names = F)


## 2. Using multinomial ENET (Only coding genes)
id_coding<-match(annotation$Ensemble_id[which(annotation$type_gene=="protein_coding")],rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_coding),]
###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog_coding),clin,family="multinomial",type.multinomial = "grouped"
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

enet<-glmnet(t(norm_data_rlog_coding),clin,family="multinomial",type.multinomial = "grouped",standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta[[2]])[which(enet$beta[[2]]!=0)]
coef1<-enet$beta[[1]][which(enet$beta[[1]]!=0)]
coef2<-enet$beta[[2]][which(enet$beta[[2]]!=0)]
results<-annotation[match(genes,annotation$Ensemble_id),]

write.csv(cbind(results,coef1,coef2),"Results/RNAseq/results_AMR_CMR_STA_ENET_coding.csv",row.names = F)



###Plot the results
results<-read.csv("Results/RNAseq/results_AMR_CMR_STA_ENET_coding.csv")
results_coding<-results[which(results$type_gene=="protein_coding"),]

id_gene<-match(results_coding$Ensemble_id,rownames(norm_data_rlog))
significantResults<-norm_data_rlog[na.omit(id_gene),]

id<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
ann_colors = list(Type = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]))

tiff(filename = "Results/RNAseq/heatmap_ENET_multinomial_coding.tiff", width = 3000, height = 2000, res = 300)
out<-pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()

cluster<-sort(cutree(out$tree_row, k=3))
id_clust<-match(results$name,names(cluster))    
results$cluster<-cluster[id_clust]
plot(out$tree_row)
abline(h=9, col="red", lty=2, lwd=2)

hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete")
cluster_similarity(hc.cut$cluster,clin) ##0.68
##clsuter plot
fviz_cluster(hc.cut, t(significantResults),ellipse.type = "norm") + theme_minimal()
# # Visualize dendrogram
##clsuter plot
tiff(filename = "Results/RNAseq/cluster_ENET_multinomial_coding.tiff", width = 2000, height = 2000, res = 300)
fviz_cluster(hc.cut, t(significantResults),ellipse.type = "norm",habillage=clin,palette=COLOR[c(2,1,3)]) + theme_minimal()
dev.off()
# # Visualize dendrogram
fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE,habillage=clin,palette=COLOR[c(2,1,3)])
# # Visualize silhouhette information
tiff(filename = "Results/RNAseq/silhouette_ENET_multinomial_coding.tiff", width = 2000, height = 2500, res = 300)
fviz_silhouette(hc.cut,palette=COLOR[c(2,1,3)])
dev.off()



