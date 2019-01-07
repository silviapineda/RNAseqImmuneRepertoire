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
coldata<-read.csv("Data/colData.csv")
coldata_rej<-read.csv("Data/colDataRej.txt")

###To study the three groups
rownames(coldata)<-coldata[,1]
coldata<-coldata[,-1]
coldata<-coldata[order(rownames(coldata)),]
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type) ##60,498

##TO study the rejectors group
rownames(coldata_rej)<-coldata_rej[,1]
coldata<-coldata_rej[,-1]
coldata<-coldata[order(rownames(coldata)),]
dds <- DESeqDataSetFromMatrix(count_rsem_genes, coldata, ~ type) ##60,498

#filtering by low counts 
##keep genes that have counts per milion values above 0.5 in at least two libraries
dds_filter <-  dds[rowSums(fpm(dds) > 0.5) >= 2] #26,545
dds_filter$type <- relevel(dds_filter$type, ref = "STA")
##normalized with rlog
normalized_rlog <- rlog(dds_filter, blind=F) 
norm_data_rlog<-assay(normalized_rlog)

##save results
clin<-ifelse(clin=="AMR", "AMR",
             ifelse(clin=="CMR","TCMR","STA"))
save(norm_data_rlog,clin,annotation,file="Data/norm_data_rlog_filter.Rdata")

clin2<-ifelse(clin=="AMR" | clin=="TCMR", "REJ","STA")
save(norm_data_rlog,clin2,annotation,file="Data/norm_data_rlog_filter_REJSTA.Rdata")

load("Data/norm_data_rlog_filter.Rdata")

###To find which ones are coding and which ones are non coding
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id_genes<-match(annotation_coding$Ensemble_id,rownames(dds_filter))
dds_filter_coding<-dds_filter[na.omit(id_genes),] #only coding
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_genes),] #15420

###Delete the Chromosome M
annotation<-annotation[which(annotation$Chr!="chrM"),]
annotation$Chr<-factor(annotation$Chr)
id<-match(rownames(norm_data_rlog_coding),annotation$Ensemble_id)
norm_data_rlog_coding<-norm_data_rlog_coding[which(is.na(id)==F),]
id<-match(rownames(dds_filter_coding),annotation$Ensemble_id)
dds_filter_coding_2<-dds_filter_coding[which(is.na(id)==F),] #15407


dds <- DESeq(dds_filter_coding_2) # N.B. Takes ~ 2 min.
res1 <- results(dds, contrast=c("type","ABMR", "STA"), alpha = .05) #length(which(res1$padj<0.05)) 4221
res2 <- results(dds, contrast=c("type","TCMR","STA"), alpha = .05) #length(which(res2$padj<0.05)) 0
res3 <- results(dds, contrast=c("type","ABMR","TCMR"), alpha = .05) #length(which(res3$padj<0.05)) 2,092

dds<-DESeq(dds_filter_coding_2)
res4 <- results(dds, contrast=c("type","REJ","STA"), alpha = .05) #length(which(res3$padj<0.05)) 2095

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
results_STA_AMR_DEseq<-annotation_sign #4221
results_STA_AMR_DEseq_FC<-results_STA_AMR_DEseq[which(abs(results_STA_AMR_DEseq$log2FC)>0.585),] #3541
write.csv(results_STA_AMR_DEseq, "Results/RNAseq/results_STA_AMR_DEseq_coding.csv", row.names = F)
write.csv(results_STA_AMR_DEseq_FC, "Results/RNAseq/results_STA_AMR_DEseq_FC_coding.csv", row.names = F)

COLOR = brewer.pal(4,"Pastel1")

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_Dseq2_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "STA", results_STA_AMR_DEseq, color=COLOR[c(1,3)])
dev.off()

#2. ENET with binomial distribution
results_STA_AMR_ENET<-ENET_binomial(clin, "AMR", "STA", norm_data_rlog_coding) ##50 genes
write.csv(results_STA_AMR_ENET,"Results/RNAseq/results_STA_AMR_ENET_coding.csv",row.names = F)

##load the results
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "STA", results_STA_AMR_ENET, color=COLOR[c(1,3)])
dev.off()

id_overlap<-match(results_STA_AMR_ENET$name,results_STA_AMR_DEseq$name)
length(which(is.na(id_overlap)==FALSE)) #45

###cluster and similarity
results_STA_AMR_ENET<-read.csv("Results/RNAseq/results_STA_AMR_ENET_coding.csv")
significantResults <- norm_data_rlog_coding[match(results_STA_AMR_ENET$Ensemble_id,rownames(norm_data_rlog_coding)),]
id_gene<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id_gene]
significantResults<-significantResults[,which(clin!="TCMR")]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
clin_2<-clin[which(clin!="TCMR")]
hc.cut <- hcut(t(xtst), k = 2, hc_method = "complete",habillage=clin_2)
library("clusteval")
clin_2<-as.numeric(data.frame(clin_2)[,1])
cluster_similarity(hc.cut$cluster,clin_2)

##################
### CMR vs STA ###
##################
results_CMR_STA_ENET<-ENET_binomial(clin, "STA", "TCMR", norm_data_rlog_coding) ##1 gene


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
results_CMR_AMR_DEseq<-annotation_sign #2092
results_CMR_AMR_DEseq_FC<-results_CMR_AMR_DEseq[which(abs(results_CMR_AMR_DEseq$log2FC)>0.585),] #1739
write.csv(results_CMR_AMR_DEseq, "Results/RNAseq/results_CMR_AMR_DEseq_coding.csv", row.names = F) 
write.csv(results_CMR_AMR_DEseq_FC, "Results/RNAseq/results_CMR_AMR_DEseq_FC_coding.csv", row.names = F)

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_CMR_AMR_Dseq2_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "TCMR", results_CMR_AMR_DEseq, color=COLOR[c(1,2)])
dev.off()

#2. ENET with binomial distribution
results_CMR_AMR_ENET<-ENET_binomial(clin, "AMR", "TCMR", norm_data_rlog_coding) ##1073 genes
write.csv(results_CMR_AMR_ENET,"Results/RNAseq/results_CMR_AMR_ENET_coding.csv",row.names = F)

tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "TCMR", results_CMR_AMR_ENET, color=COLOR[c(1,2)])
dev.off()

id_overlap<-match(results_CMR_AMR_ENET$name,results_CMR_AMR_DEseq$name)
length(which(is.na(id_overlap)==FALSE)) #538

###cluster and similarity
results_CMR_AMR_ENET<-read.csv("Results/RNAseq/results_CMR_AMR_DEseq_FC_coding.csv")
significantResults <- norm_data_rlog_coding[match(results_CMR_AMR_ENET$Ensemble_id,rownames(norm_data_rlog_coding)),]
id_gene<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id_gene]
significantResults<-significantResults[,which(clin!="STA")]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
clin_2<-clin[which(clin!="STA")]
hc.cut <- hcut(t(xtst), k = 2, hc_method = "complete",habillage=clin_2)
library("clusteval")
clin_2<-as.numeric(data.frame(clin_2)[,1])
cluster_similarity(hc.cut$cluster,clin_2)

##################
### Rej vs. STA ###
##################
#results
allgenes <- rownames(res4)
genes_sign <- allgenes[which(res4$padj < 0.05)]
id_sign<-match(genes_sign,annotation$Ensemble_id)
annotation_sign<-annotation[id_sign,] 
annotation_sign$log2FC<-res4$log2FoldChange[which(res4$padj < 0.05)]
results_STA_REJ_DEseq<-annotation_sign #1094
results_STA_REJ_DEseq_FC<-results_STA_REJ_DEseq[which(abs(results_STA_REJ_DEseq$log2FC)>0.585),] #875
write.csv(results_STA_REJ_DEseq, "Results/RNAseq/results_STA_REJ_DEseq_coding.csv", row.names = F)
write.csv(results_STA_REJ_DEseq_FC, "Results/RNAseq/results_STA_REJ_DEseq_coding_FC.csv", row.names = F)
COLOR = brewer.pal(4,"Pastel1")

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_REJ_STA_Dseq2_coding_FC.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin2,norm_data_rlog, "REJ", "STA", results_STA_REJ_DEseq_FC, color=COLOR[c(4,1)])
dev.off()

#2. ENET with binomial distribution
results_REJ_STA_ENET<-ENET_binomial(clin2, "REJ", "STA", norm_data_rlog_coding) ##328 genes
write.csv(results_REJ_STA_ENET,"Results/RNAseq/results_REJ_STA_ENET_coding.csv",row.names = F)

tiff(filename = "Results/RNAseq/heatmap_REJ_STA_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
out<-plot_heatmap(clin2,norm_data_rlog, "REJ", "STA", results_REJ_STA_ENET, color=COLOR[c(4,1)])
dev.off()

id_overlap<-match(results_REJ_STA_ENET$name,results_STA_REJ_DEseq$name)
length(which(is.na(id_overlap)==FALSE)) #205

###cluster and similarity
results_REJ_STA_ENET<-read.csv("Results/RNAseq/results_REJ_STA_ENET_coding.csv")
significantResults <- norm_data_rlog_coding[match(results_REJ_STA_ENET$Ensemble_id,rownames(norm_data_rlog_coding)),]
id_gene<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id_gene]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)

hc.cut <- hcut(t(xtst), k = 2, hc_method = "complete",habillage=clin2)
library("clusteval")
clin3<-as.numeric(data.frame(clin2)[,1])
cluster_similarity(hc.cut$cluster,clin3)


#########################
## AMR vs CMR vs STA ###
########################

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

write.csv(cbind(results,coef1,coef2),"Results/RNAseq/results_AMR_CMR_STA_ENET_coding.csv",row.names = F) #131 

###Plot the results
results<-read.csv("Results/RNAseq/results_AMR_CMR_STA_ENET_coding.csv")

id_gene<-match(results$Ensemble_id,rownames(norm_data_rlog_coding))
significantResults<-norm_data_rlog_coding[na.omit(id_gene),]

id<-match(rownames(significantResults),annotation$Ensemble_id)
rownames(significantResults)<-annotation$name[id]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Pastel1")
ann_colors = list(Type = c("AMR" = COLOR [1] ,"TCMR" = COLOR[2], "STA" = COLOR[3]))

tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_STA_ENET_coding.tiff", width = 3000, height = 3500, res = 300)
out<-pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()

###Clusters
cluster<-sort(cutree(out$tree_row, k=3))
id_clust<-match(results$name,names(cluster))    
results$cluster<-cluster[id_clust]
write.csv(results,"Results/RNAseq/results_AMR_CMR_STA_ENET_coding.csv",row.names = F)  
plot(out$tree_row)
abline(h=9, col="red", lty=2, lwd=2)



## Obtain the clusterization
#COLOR = add.alpha(brewer.pal(4,"Set1"),0.3)
COLOR = brewer.pal(4,"Pastel1")
hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete",habillage=clin)
clin_3<-as.numeric(data.frame(clin)[,1])
cluster_similarity(hc.cut$cluster,clin_3) ##1
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


######Gene Ontology
library(AnnotationHub)
hub <- AnnotationHub()
