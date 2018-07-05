rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: QC analysis
###           
### DESCRIP: Analysis of the output from the alignment
###         
###
### Author: Silvia Pineda
### Date: October, 2017
############################################################################################
library(glmnet)
library("pheatmap")
library("plyr")
library("RColorBrewer")

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/repertoireResults.Rdata")
load("Data/RNAseqProcessedNormalized.Rdata")
load("Data/norm_data_rlog_filter.Rdata")

###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog),summaryMatrix$Alpha_Beta_ratio_expression,standardize=TRUE,
                                                        alpha=a,nfolds=5))})
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

enet<-glmnet(t(norm_data_rlog),summaryMatrix$Alpha_Beta_ratio_expression,standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]
results<-annotation[match(genes,annotation$Ensemble_id),]
write.csv(cbind(results,coef),"Results/Integration/genes.enet.alphaOptimum.csv",row.names = F)

###############
## Heatmap ####
##############
significant_norm<-norm_data_rlog[match(genes,rownames(norm_data_rlog)),]
rownames(significant_norm)<-results$name

xt<-t(significant_norm)
xts<-scale(xt)
xtst<-t(xts)

annotation_col = data.frame(
  phenotype = clin,
  ratio = summaryMatrix$Alpha_Beta_ratio_expression)

rownames(annotation_col)<-colnames(significant_norm)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]),
                   ratio = brewer.pal(6,"Reds"))
tiff("Results/Integration/heatmap_ratio_integration.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_rows = T, cluster_cols = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,border_color=F)
dev.off()


###Restricted to protein_coding
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id.genes<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id.genes),] ## genes annotated 15,420

fit<-VSURF(x = t(norm_data_rlog_coding), y=summaryMatrix$Alpha_Beta_ratio_expression, parallel = TRUE,ncores=4)
save(fit,file="Results/RNAseq/Genes_VSURF_alphaBeta.Rdata")


alphalist<-seq(0.01,0.99,by=0.01)
set.seed(33)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog_coding),summaryMatrix$Alpha_Beta_ratio_expression,standardize=TRUE,
                                                        alpha=a,nfolds=5))})
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

enet<-glmnet(t(norm_data_rlog_coding),summaryMatrix$Alpha_Beta_ratio_expression,standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta)[which(as.numeric(enet$beta)!=0)]
coef<-enet$beta[which(as.numeric(enet$beta)!=0)]
results<-cbind(annotation_coding[match(genes,annotation_coding$Ensemble_id),],coef)
write.csv(cbind(results,coef),"Results/Integration/genes.enet.alphaOptimum.coding.csv",row.names = F)



###############
## Heatmap ####
##############
significant_norm_coding<-norm_data_rlog_coding[match(genes,rownames(norm_data_rlog_coding)),]
rownames(significant_norm_coding)<-results$name

xt<-t(significant_norm_coding)
xts<-scale(xt)
xtst<-t(xts)

annotation_col = data.frame(
  phenotype = clin,
  ratio = summaryMatrix$Alpha_Beta_ratio_expression)

rownames(annotation_col)<-colnames(significant_norm_coding)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]),
                   ratio = brewer.pal(6,"Reds"))
tiff("Results/Integration/heatmap_ratio_coding_integration.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_cols = T,cluster_rows = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,border_color=F)
dev.off()


###Study genes
COLOR=brewer.pal(3,"Set2")
id<-grep("PPCD",rownames(significant_norm_coding))
plot(significant_norm_coding[id,],summaryMatrix$Alpha_Beta_ratio_expression,col=COLOR[clin],pch=19)
boxplot(significant_norm_coding[id,]~clin,col=COLOR)


##########################
## Run RF in the server ##
##########################
load("repertoireResults.Rdata")
load("norm_data_rlog_filter.Rdata")

annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id.genes<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id.genes),] ## genes annotated 15,420
set.seed(33)
library(VSURF)
fit<-VSURF(x = t(norm_data_rlog_coding), y=summaryMatrix$Alpha_Beta_ratio_expression, parallel = TRUE,ncores=10)
save(fit,file="Genes_VSURF_alphaBeta.Rdata")


load("Results/Integration/Genes_VSURF_alphaBeta_NonCoding.Rdata")
genes<-data.frame(norm_data_rlog[fit$varselect.interp,]) ##8 genes

id<-match(rownames(genes),annotation$Ensemble_id)
rownames(genes)<-annotation$name[id]

xt<-t(genes)
xts<-scale(xt)
xtst<-t(xts)

annotation_col = data.frame(
  phenotype = clin,
  ratio = summaryMatrix$Alpha_Beta_ratio_expression)

rownames(annotation_col)<-colnames(genes)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]),
                   ratio = brewer.pal(6,"Reds"))
tiff("Results/Integration/heatmap_ratio_NonCoding_RF.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_cols = T,cluster_rows = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,border_color=F)
dev.off()


###Study genes
COLOR=brewer.pal(3,"Set2")
id<-grep("PPCD",rownames(genes))
plot(as.numeric(as.character(genes[id,])),summaryMatrix$Alpha_Beta_ratio_expression,col=COLOR[clin],pch=19)
boxplot(as.numeric(as.character(genes[id,]))~clin,col=COLOR)

