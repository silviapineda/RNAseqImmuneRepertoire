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
#######################
### Reproducibility ##
######################
results<-read.csv("Results/RNAseq/results_AMR_CMR_STA_ENET.csv")
genes<-list()
for(i in 35:100){
  print(i)
  alphalist<-seq(0.01,0.99,by=0.01)
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
  genes_ensembl<-rownames(enet$beta[[2]])[which(enet$beta[[2]]!=0)]
  genes[[i]]<-as.character(annotation[match(genes_ensembl,annotation$Ensemble_id),"name"])
}
