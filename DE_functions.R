rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### Functions for the DE_analysis
###         
###
### Author: Silvia Pineda
### Date: June, 2018
############################################################################################


ENET_binomial <- function(clin, categ1, categ2, norm_data_rlog){
 
  clin_2<-factor(clin[which(clin==categ1 | clin==categ2)])
  norm_data_rlog_2<-norm_data_rlog[,which(clin==categ1 | clin==categ2)]
 
  ##CV to find the best alpha and lambda 
  alphalist<-seq(0.01,0.99,by=0.01)
  set.seed(54)
  elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog_2),clin_2,family="binomial",
                                                        standardize=TRUE,alpha=a,nfolds=5))})
  xx<-rep(NA,length(alphalist))
  yy<-rep(NA,length(alphalist))
  for (j in 1:length(alphalist)) {
    if(class(elasticnet[[j]]) != "try-error"){
      xx[j]<-elasticnet[[j]]$lambda.min
      id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
      yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
    }
  }
  id.min<-which(yy==min(yy,na.rm=TRUE))
  lambda<-xx[id.min]
  alpha<-alphalist[id.min]

  enet<-glmnet(t(norm_data_rlog_2),clin_2,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda)
  genes<-rownames(enet$beta)[which(enet$beta!=0)]
  coef<-enet$beta[which(enet$beta!=0)]
  results<-annotation[match(genes,annotation$Ensemble_id),]
  
  ###Obtain the FC
  id_sign<-match(results$Ensemble_id,rownames(norm_data_rlog_2))
  if(length(id_sign)>1){
    results$mean1<-apply(norm_data_rlog_2[id_sign,which(clin_2==categ1)],1,mean)
    results$mean2<-apply(norm_data_rlog_2[id_sign,which(clin_2==categ2)],1,mean)
  } else {
    results$mean1<-mean(norm_data_rlog_2[id_sign,which(clin_2==categ1)])
    results$mean2<-mean(norm_data_rlog_2[id_sign,which(clin_2==categ2)])
  }
  results$Log2FC<-results$mean1-results$mean2
  return(cbind(results,coef))
} 

LASSO_binomial <- function(clin, categ1, categ2, norm_data_rlog){
  
  clin_2<-factor(clin[which(clin==categ1 | clin==categ2)])
  norm_data_rlog_2<-norm_data_rlog[,which(clin==categ1 | clin==categ2)]
  
  set.seed(54)
  lassocv<-cv.glmnet(t(norm_data_rlog_2),clin_2,family="binomial",standardize=TRUE,alpha=1,nfolds=5)
  lasso<-glmnet(t(norm_data_rlog_2),clin_2,family="binomial",standardize=TRUE,alpha=1,lambda=lassocv$lambda)
  genes<-rownames(enet$beta)[which(enet$beta!=0)]
  coef<-enet$beta[which(enet$beta!=0)]
  results<-annotation[match(genes,annotation$Ensemble_id),]
  return(cbind(results,coef))
} 


plot_heatmap <- function(clin,norm_data_rlog, categ1, categ2, results, color){
  
  clin_2<-factor(clin[which(clin==categ1 | clin==categ2)])
  norm_data_rlog_2<-norm_data_rlog[,which(clin==categ1 | clin==categ2)]
  significantResults <- norm_data_rlog_2[match(results$Ensemble_id,rownames(norm_data_rlog_2)),]
  id_gene<-match(rownames(significantResults),annotation$Ensemble_id)
  rownames(significantResults)<-annotation$name[id_gene]

  xt<-t(significantResults)
  xts<-scale(xt)
  xtst<-t(xts)
  rownames <- colnames(significantResults)
  annotation.col <- data.frame(row.names = rownames)
  annotation.col$Type <- factor(clin_2)
  ann_colors = list(Type = color)
  names(ann_colors$Type)<-levels(clin_2)
  p<-pheatmap(xtst, annotation = annotation.col, annotation_colors = ann_colors,show_rownames=T,border_color=F)
  return(p)
}



