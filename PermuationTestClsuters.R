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
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("glmnet")
library("GISTools")
library(factoextra)
require("cluster")

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

load("Data/RNAseqProcessedNormalized.Rdata")
load("Data/norm_data_rlog_filter.Rdata")
###############################################
### permutation analysis for the clustering ##
##############################################

## 2. Using multinomial ENET
set.seed(54)
similarity<-NULL
for (i in 1:10){
  print(i)
  rm(clin_perm)
  clin_perm<-sample(clin)
  write.csv(clin_perm,paste0("Results/RNAseq/clin_perm",i,".txt"))
  
  ###1. Apply ENET
  alphalist<-seq(0.01,0.99,by=0.01)
  elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog),clin_perm,family="multinomial",type.multinomial = "grouped"
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
  
  enet<-glmnet(t(norm_data_rlog),clin_perm,family="multinomial",type.multinomial = "grouped",standardize=TRUE,alpha=alpha,lambda=lambda)
  genes<-rownames(enet$beta[[2]])[which(enet$beta[[2]]!=0)]
  coef1<-enet$beta[[1]][which(enet$beta[[1]]!=0)]
  coef2<-enet$beta[[2]][which(enet$beta[[2]]!=0)]
  results<-annotation[match(genes,annotation$Ensemble_id),]
  write.csv(results,paste0("Results/RNAseq/perm",i,".csv"),row.names = F)
  
  ##2. plot heatmap
  if(dim(results)[1]>1){
    id_gene<-match(results$Ensemble_id,rownames(norm_data_rlog))
    significantResults<-norm_data_rlog[na.omit(id_gene),]
    
    id<-match(rownames(significantResults),annotation$Ensemble_id)
    rownames(significantResults)<-annotation$name[id]
    
    xt<-t(significantResults)
    xts<-scale(xt)
    xtst<-t(xts)
    rownames <- colnames(significantResults)
    annotation.col <- data.frame(row.names = rownames)
    annotation.col$Type <- factor(clin_perm)
    COLOR = brewer.pal(4,"Pastel1")
    ann_colors = list(Type = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]))
    
    tiff(filename = paste0("Results/RNAseq/heatmap_perm",i,".tiff"), width = 3000, height = 2000, res = 300)
    pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
    dev.off()
    
    #3. Obtain similarity measure
    hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete")
    similarity[i]<-cluster_similarity(hc.cut$cluster,clin_perm) 
  } else {
    similarity[i]<-0
  }
}
mean(similarity) ##0.27


obtain_cluster_perm <- function(xtst,significantResults){
  COLOR = brewer.pal(4,"Pastel1")
  ## Obtain the clusterization
  hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete")
  ##clsuter plot
  fviz_cluster(hc.cut, t(significantResults),ellipse.type = "norm",palette=COLOR[c(2,1,3)]) + theme_minimal()
  return(hc.cut$cluster)
}

obtain_sil_perm <- function(xtst){
  COLOR = brewer.pal(4,"Pastel1")
  ## Obtain the clusterization
  hc.cut <- hcut(t(xtst), k = 3, hc_method = "complete")
  fviz_silhouette(hc.cut,palette=COLOR[c(2,1,3)])
}

