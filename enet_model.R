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
library("VSURF")
library("randomForest")

setwd("~/ImmuneRep_RNAseq/Integration/")
load("/Users/Pinedasans/ImmuneRep_RNAseq/MIXCR/repertoireResults.Rdata")
load("/Users/Pinedasans/ImmuneRep_RNAseq/QC/rsemResults.Rdata")
load("/Users/Pinedasans/ImmuneRep_RNAseq/Integration/norm_data_filter_rlog.Rdata")

annotation<-read.table("/Users/Pinedasans/ImmuneRep_RNAseq/Integration/gencode.v23.annotation.txt",sep=";") ##60,498
annotation<-annotation[,-7]
colnames(annotation)<-c("Chr","Start","End","Ensemble_id","type_gene","name")

###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(1234)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog),repertoireResults$Alpha_Beta,standardize=TRUE,
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

enet<-glmnet(t(norm_data_rlog),repertoireResults$Alpha_Beta,standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]
results<-annotation[match(genes,annotation$Ensemble_id),]
write.csv(cbind(results,coef),"genes.enet.alphaOptimum.csv",row.names = F)

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
  ratio = repertoireResults$Alpha_Beta)

rownames(annotation_col)<-colnames(significant_norm)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]),
                   ratio = brewer.pal(6,"Reds"))
tiff("heatmap_ratio_noncoding.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_rows = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,border_color=F)
dev.off()


###Restricted to protein_coding
annotation_coding<-annotation[which(annotation$type_gene=="protein_coding"),]
id.genes<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id.genes),] ## genes annotated 17,832

alphalist<-seq(0.01,0.99,by=0.01)
set.seed(99999)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog_coding),repertoireResults$Alpha_Beta,standardize=TRUE,
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

enet<-glmnet(t(norm_data_rlog_coding),repertoireResults$Alpha_Beta,standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta)[which(as.numeric(enet$beta)!=0)]
coef<-enet$beta[which(as.numeric(enet$beta)!=0)]
results<-cbind(annotation_coding[match(genes,annotation_coding$Ensemble_id),],coef)
write.csv(cbind(results,coef),"genes.enet.alphaOptimum.coding.csv",row.names = F)



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
  ratio = repertoireResults$Alpha_Beta)

rownames(annotation_col)<-colnames(significant_norm_coding)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]),
                   ratio = brewer.pal(6,"Reds"))
tiff("heatmap_ratio_coding_opt.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_cols = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,border_color=F)
dev.off()


#####################
### Random Forest ###
#####################
set.seed(112233)

fit<-VSURF(x = t(norm_data_rlog_coding), y=clin,parallel = TRUE,ncores=4)
save(fit,file="CodingGenes_RF.Rdata")
genes<-data.frame(norm_data_rlog_coding[fit$varselect.interp,])

set.seed(1000)
rf_output <- randomForest(clin~.,data=t(genes),proximity=TRUE, keep.forest=T,ntree=1000)

##Heatrmap
significant_norm_coding<-genes

xt<-t(significant_norm_coding)
xts<-scale(xt)
xtst<-t(xts)

annotation_col = data.frame(
  phenotype = clin)

rownames(annotation_col)<-colnames(significant_norm_coding)
fill=brewer.pal(3,"Set2")
ann_colors = list (phenotype = c(STA = fill[1], CMR = fill[2], AMR = fill[3]))
tiff("heatmap_DE_coding_RF.tiff",res=300,w=3000,h=2000)
pheatmap(xtst,cluster_cols = T,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,border_color=F)
dev.off()
