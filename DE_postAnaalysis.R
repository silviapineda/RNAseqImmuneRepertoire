rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Post DE analysis
###           
### DESCRIP: Analysis with the data from the QC process
###         
###
### Author: Silvia Pineda
### Date: October, 2017
############################################################################################
library("RColorBrewer")
library("igraph")
library("qgraph")
library("corrplot")
library("gplots")
library("pheatmap")
library(dplyr)

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

###Load the data
load("Data/norm_data_rlog_filter.Rdata")

results_multi<-read.csv("Results/RNAseq/results_AMR_CMR_STA_ENET.csv")
results_AMR_STA<-read.csv("Results/RNAseq/results_STA_AMR_ENET.csv")
results_AMR_CMR<-read.csv("Results/RNAseq/results_CMR_AMR_ENET.csv")
results_REJ_STA<-read.csv("Results/RNAseq/results_REJ_STA_ENE_codingT.csv")
results_REJ_STA_DESeq<-read.csv("Results/RNAseq/results_STA_REJ_DEseq.csv")

##### Clinical Data Analysis #######
clinData<-read.csv("Data/ClinicalData.csv")
clinData<-clinData[order(clinData$Individual_id),]
clinData<-clinData[which(clinData$Individual_id!="M1" & clinData$Individual_id!="M2" & clinData$Individual_id!="M3" & clinData$Individual_id!="M4"),]
clinData$clin<-c(rep("CMR",13),rep("AMR",12),rep("STA",12))

##### Demographics ####
#donor age
mean(clinData$donor.age[which(clinData$clin=="CMR")])
sd(clinData$donor.age[which(clinData$clin=="CMR")])
mean(clinData$donor.age[which(clinData$clin=="AMR")])
sd(clinData$donor.age[which(clinData$clin=="AMR")])
mean(clinData$donor.age[which(clinData$clin=="STA")])
sd(clinData$donor.age[which(clinData$clin=="STA")])
summary(lm(clinData$donor.age~clinData$clin))

#recipient age
mean(clinData$Rec.age[which(clinData$clin=="CMR")])
sd(clinData$Rec.age[which(clinData$clin=="CMR")])
mean(clinData$Rec.age[which(clinData$clin=="AMR")])
sd(clinData$Rec.age[which(clinData$clin=="AMR")])
mean(clinData$Rec.age[which(clinData$clin=="STA")])
sd(clinData$Rec.age[which(clinData$clin=="STA")])
summary(lm(clinData$Rec.age~clinData$clin))

#recipient gender
table(clinData$Rec.Gender,clinData$clin)
chisq.test(table(clinData$Rec.Gender,clinData$clin))

#donor gender
table(clinData$Donor..gender,clinData$clin)
chisq.test(table(clinData$Donor..gender,clinData$clin))

#Number of TX
mean(clinData$number.TX[which(clinData$clin=="CMR")])
sd(clinData$number.TX[which(clinData$clin=="CMR")])
mean(clinData$number.TX[which(clinData$clin=="AMR")])
sd(clinData$number.TX[which(clinData$clin=="AMR")])
mean(clinData$number.TX[which(clinData$clin=="STA")])
sd(clinData$number.TX[which(clinData$clin=="STA")])
summary(lm(clinData$number.TX~clinData$clin))

#Induction type
table(clinData$Induction..Bxb..1..or.rATg..2..none..0..,clinData$clin)
chisq.test(table(clinData$Induction..Bxb..1..or.rATg..2..none..0..,clinData$clin))

chisq.test(table(clinData$CNI.time.rejection..FK.1.CsA.2.none.0.,clinData$clin))

#DSA
table(clinData$DSA_HLA.class..0.neg..1.class.I..2.class.II..3.class.I.and.II.,clinData$clin)
chisq.test(table(clinData$DSA_HLA.class..0.neg..1.class.I..2.class.II..3.class.I.and.II.,clinData$clin))

##Anti-HLA
table(clinData$Anti.HLA.Ab_HLA.class,clinData$clin)
chisq.test(table(clinData$Anti.HLA.Ab_HLA.class,clinData$clin))

#eGFR
mean(clinData$eGFR..time.Bx.[which(clinData$clin=="CMR")])
sd(clinData$eGFR..time.Bx.[which(clinData$clin=="CMR")])
mean(clinData$eGFR..time.Bx.[which(clinData$clin=="AMR")])
sd(clinData$eGFR..time.Bx.[which(clinData$clin=="AMR")])
mean(clinData$eGFR..time.Bx.[which(clinData$clin=="STA")])
sd(clinData$eGFR..time.Bx.[which(clinData$clin=="STA")])
summary(lm(clinData$eGFR..time.Bx.~clinData$clin))

#GraftLoss
table(clinData$GraftLossCateg,clinData$clin)
chisq.test(table(clinData$GraftLossCateg,clinData$clin))

#Immunosuppression
table(clinData$Immunosuppression,clinData$clin)

summary(glm(summaryMatrix$Alpha_Beta_ratio_expression~clinData$clin+clinData$Immunosuppression)) #
boxplot(summaryMatrix$Alpha_Beta~clinData$donor.age,col=COLOR)

#ag
mean(clinData$ag..glomerulitis.[which(clinData$clin=="CMR")])
sd(clinData$ag..glomerulitis.[which(clinData$clin=="CMR")])
mean(clinData$ag..glomerulitis.[which(clinData$clin=="AMR")])
sd(clinData$ag..glomerulitis.[which(clinData$clin=="AMR")])
mean(clinData$ag..glomerulitis.[which(clinData$clin=="STA")])
sd(clinData$ag..glomerulitis.[which(clinData$clin=="STA")])
summary(lm(clinData$ag..glomerulitis.~clinData$clin))
#ai
mean(clinData$ai..interstitial.inflammation.[which(clinData$clin=="CMR")])
sd(clinData$ai..interstitial.inflammation.[which(clinData$clin=="CMR")])
mean(clinData$ai..interstitial.inflammation.[which(clinData$clin=="AMR")])
sd(clinData$ai..interstitial.inflammation.[which(clinData$clin=="AMR")])
mean(clinData$ai..interstitial.inflammation.[which(clinData$clin=="STA")])
sd(clinData$ai..interstitial.inflammation.[which(clinData$clin=="STA")])
summary(lm(clinData$ai..interstitial.inflammation.~clinData$clin))
#at
mean(clinData$at..tubulitis.[which(clinData$clin=="CMR")])
sd(clinData$at..tubulitis.[which(clinData$clin=="CMR")])
mean(clinData$at..tubulitis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$at..tubulitis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$at..tubulitis.[which(clinData$clin=="STA")])
sd(clinData$at..tubulitis.[which(clinData$clin=="STA")])
summary(lm(clinData$at..tubulitis.~clinData$clin))
#ti
mean(clinData$ti..total.inflammation.[which(clinData$clin=="CMR")])
sd(clinData$ti..total.inflammation.[which(clinData$clin=="CMR")])
mean(clinData$ti..total.inflammation.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$ti..total.inflammation.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$ti..total.inflammation.[which(clinData$clin=="STA")])
sd(clinData$ti..total.inflammation.[which(clinData$clin=="STA")])
summary(lm(clinData$ti..total.inflammation.~clinData$clin))

#ptc
mean(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="CMR")])
sd(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="CMR")])
mean(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="STA")])
sd(clinData$ptc..peritubular.capilaritis.[which(clinData$clin=="STA")])
summary(lm(clinData$ptc..peritubular.capilaritis.~clinData$clin))
#av
mean(clinData$av..endothelialitis.[which(clinData$clin=="CMR")])
sd(clinData$av..endothelialitis.[which(clinData$clin=="CMR")])
mean(clinData$av..endothelialitis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$av..endothelialitis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$av..endothelialitis.[which(clinData$clin=="STA")])
sd(clinData$av..endothelialitis.[which(clinData$clin=="STA")])
summary(lm(clinData$av..endothelialitis.~clinData$clin))
#C4D
mean(clinData$C4d[which(clinData$clin=="CMR")])
sd(clinData$C4d[which(clinData$clin=="CMR")])
mean(clinData$C4d[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$C4d[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$C4d[which(clinData$clin=="STA")])
sd(clinData$C4d[which(clinData$clin=="STA")])
summary(lm(clinData$C4d~clinData$clin))

##Chronic
#cg
mean(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="CMR")])
sd(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="CMR")])
mean(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="STA")])
sd(clinData$cg..transplant.glomerulopathy.[which(clinData$clin=="STA")])
summary(lm(clinData$cg..transplant.glomerulopathy.~clinData$clin))
#ci
mean(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="CMR")])
sd(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="CMR")])
mean(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="STA")])
sd(clinData$ci..intestitial.fibrosis.[which(clinData$clin=="STA")])
summary(lm(clinData$ci..intestitial.fibrosis.~clinData$clin))
#ct
mean(clinData$ct..tubular.atrophy.[which(clinData$clin=="CMR")])
sd(clinData$ct..tubular.atrophy.[which(clinData$clin=="CMR")])
mean(clinData$ct..tubular.atrophy.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$ct..tubular.atrophy.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$ct..tubular.atrophy.[which(clinData$clin=="STA")])
sd(clinData$ct..tubular.atrophy.[which(clinData$clin=="STA")])
summary(lm(clinData$ct..tubular.atrophy.~clinData$clin))
#cv
mean(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="CMR")])
sd(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="CMR")])
mean(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="STA")])
sd(clinData$cv.endothelial.fibrosis.[which(clinData$clin=="STA")])
summary(lm(clinData$cv.endothelial.fibrosis.~clinData$clin))
#ah
mean(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="CMR")])
sd(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="CMR")])
mean(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="STA")])
sd(clinData$aah..arteriolar.hyalinosis.[which(clinData$clin=="STA")])
summary(lm(clinData$aah..arteriolar.hyalinosis.~clinData$clin))
#cm
mean(clinData$cm..mesangial.expansion.[which(clinData$clin=="CMR")])
sd(clinData$cm..mesangial.expansion.[which(clinData$clin=="CMR")])
mean(clinData$cm..mesangial.expansion.[which(clinData$clin=="AMR")],na.rm=T)
sd(clinData$cm..mesangial.expansion.[which(clinData$clin=="AMR")],na.rm = T)
mean(clinData$cm..mesangial.expansion.[which(clinData$clin=="STA")])
sd(clinData$cm..mesangial.expansion.[which(clinData$clin=="STA")])
summary(lm(clinData$cm..mesangial.expansion.~clinData$clin))



##########################
### Enrichment Analysis ##
##########################
##AMR vs. STA
table(results_AMR_STA$type_gene)
table(results_AMR_STA$cluster)
table(results_AMR_STA$type_gene,results_AMR_STA$cluster)
### AMR in non-coding is enriched
counts=matrix(data=c(13,30,15,45),nrow=2)
fisher.test(counts) #p-value = 0.6

##Is enriched respect the selected in the STA?
results_AMR_STA$type_gene_2<-ifelse(results_AMR_STA$type_gene!="protein_coding","non-coding","coding")
tab<-table(results_AMR_STA$cluster,results_AMR_STA$type_gene_2)
fisher.test(tab) #p-value = 0.2

##AMR vs. CMR
table(results_AMR_CMR$type_gene)
table(results_AMR_CMR$cluster)
table(results_AMR_CMR$type_gene,results_AMR_CMR$cluster)
### AMR in non-coding
counts=matrix(data=c(7,8,8,13),nrow=2)
fisher.test(counts) #p-value = 0.7

##Is enriched respect the selected in the STA?
results_AMR_CMR$type_gene_2<-ifelse(results_AMR_CMR$type_gene!="protein_coding","non-coding","coding")
tab<-table(results_AMR_CMR$cluster,results_AMR_CMR$type_gene_2)
fisher.test(tab) #p-value = 0.3


##AMR vs. CMR vs. STA
table(results_multi$type_gene)
table(results_multi$cluster)
table(results_multi$type_gene,results_multi$cluster)
results_multi$type_gene_2<-ifelse(results_multi$type_gene!="protein_coding","non-coding","coding")

### AMR in non-coding
counts=matrix(data=c(16,32,39,63),nrow=2)
fisher.test(counts) #p-value = 0.8

### CMR in non-coding
counts=matrix(data=c(6,16,39,63),nrow=2)
fisher.test(counts) #p-value = 0.5

### STA in non-coding
counts=matrix(data=c(15,15,39,63),nrow=2)
fisher.test(counts) #p-value = 0.3

##Is enriched respect the selected in the STA?
tab<-table(results_multi$cluster,results_multi$type_gene_2)
fisher.test(tab) #p-value = 0.2

################################
### Overlap between results ###
###############################
results_AMR_CMR[which(results_AMR_CMR$name %in% results_AMR_STA$name),] #3 genes are in common
results_AMR_STA[which(results_AMR_STA$name %in% results_AMR_CMR$name),]

results_multi[which(results_AMR_STA$name %in% results_multi$name),]
results_multi[which(results_AMR_CMR$name %in% results_multi$name),]


###############################################
### Network analysis with clinical variables ##
###############################################
id_genes<-match(results_multi$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_results<-norm_data_rlog[id_genes,]

#put the names of the genes
id_genes<-match(rownames(norm_data_rlog_results),results_multi$Ensemble_id)
rownames(norm_data_rlog_results)<-results_multi$name[id_genes]

clin_vars<-clinData[,c(3,5:17,41)]
colnames(clin_vars)<-c("DSA","ag","ai","at","ti","ptc","av","ah","cg","ci","ct","cv","cm","C4d","Time")
rownames(clin_vars)<-clinData$Individual_id

####################################################
##  Correlation matrix between clinical variables ##

M <- cor(na.omit(clin_vars))
p.mat <- cor.mtest(na.omit(clin_vars))$p
ord <- corrMatOrder(M, order = "AOE")
M2 <- M[ord,ord]
p.mat2 <-p.mat[ord,ord]
tiff("Results/RNAseq/Corrplot_clinVars.tiff", res=300,w=2000,h=2000)
corrplot.mixed(M2, lower.col = brewer.pal(n = 8, name = "RdBu"), upper.col = brewer.pal(n = 8, name = "RdBu"))
dev.off()

tiff("Results/RNAseq/Corrplot_clinVars.tiff", res=300,w=2000,h=2000)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M2, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "red", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat2, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()

##########################
#### Linear regression ##
#########################
##1. Find all the p-values
pvalue<-matrix(NA,nrow(norm_data_rlog_results),ncol(clin_vars))
for (i in 1:nrow(norm_data_rlog_results)){
  for (j in 1:ncol(clin_vars)){
    pvalue[i,j]<-coef(summary(lm(norm_data_rlog_results[i,]~clin_vars[,j], data = clin_vars)))[2,4]
  }
}
pvalue_adj<-p.adjust(pvalue,method = "fdr")  
maxp<-max(pvalue_adj[which(pvalue_adj<0.05)])
pos<-grep(maxp, pvalue_adj) 
pvalue_array<-as.numeric(pvalue)
fdr_value<-max(pvalue_array[pos]) ##This is the FDR correction < 0.05

##2.Build the network (fdr<0.05)
lm_pvalue<-rep(NA,ncol(clin_vars))
lm_coef<-rep(NA,ncol(clin_vars))
lm_network_p<-list()
lm_network_c<-list()
for (i in 1:nrow(norm_data_rlog_results)){
  for (j in 1:ncol(clin_vars)){
    lm_pvalue[j]<-coef(summary(lm(norm_data_rlog_results[i,]~clin_vars[,j], data = clin_vars)))[2,4]
    lm_coef[j]<-coef(summary(lm(norm_data_rlog_results[i,]~clin_vars[,j], data = clin_vars)))[2,1]
  }
  sign_pvalue<-lm_pvalue[which(lm_pvalue<fdr_value)]
  sign_coef<-lm_coef[which(lm_pvalue<fdr_value)]
  names(sign_pvalue)<-colnames(clin_vars)[which(lm_pvalue<fdr_value)]
  names(sign_coef)<-colnames(clin_vars)[which(lm_pvalue<fdr_value)]
  lm_network_p[[i]]<-sign_pvalue
  lm_network_c[[i]]<-sign_coef
}
names(lm_network_p)<-rownames(norm_data_rlog_results)
names(lm_network_c)<-rownames(norm_data_rlog_results)

##Build the network
##For lm
edges<- data.frame(
  Genes = rep(names(lm_network_p), lapply(lm_network_p, length)), Clin = unlist(lapply(lm_network_p, names)),
  pvalue = abs(log(unlist(lm_network_p),10)), coef = unlist(lm_network_c))
edges$color<-ifelse(edges$coef<0,"darkgoldenrod3","gray56")

#Unique Genes with edges
Genes<-as.character(unique(edges$Genes))
#All genes significant in the DE analysis
#Genes<-as.character(results$name)
id_AMR<-na.omit(match(results_multi$name[which(results_multi$cluster=="AMR")],Genes))
id_CMR<-na.omit(match(results_multi$name[which(results_multi$cluster=="CMR")],Genes))
id_STA<-na.omit(match(results_multi$name[which(results_multi$cluster=="STA")],Genes))

nodes<-c(Genes[id_AMR],Genes[id_CMR],Genes[id_STA],colnames(clin_vars))
COLOR = brewer.pal(4,"Pastel1")
nodes<-cbind(nodes,c(rep(COLOR[1],length(id_AMR)),rep(COLOR[2],length(id_CMR)),
                     rep(COLOR[3],length(id_STA)),rep(COLOR[4],ncol(clin_vars))))
nodes<-data.frame(cbind(nodes,c(rep(5,length(Genes)),rep(10,ncol(clin_vars)))))
colnames(nodes)<-c("names","color","size")
nodes$size<-as.numeric(as.character(nodes$size))
nodes<-nodes[-90,]
#Make the graph and plot
net <- graph_from_data_frame(d = edges, vertices = nodes,directed = F)

E(net)$width<-E(net)$pvalue
tiff(paste("Results/RNAseq/network_lm_fdr_05.tiff",sep=""),res=300,h=5000,w=5000)
l <- layout_with_fr(net)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net,vertex.label.color="black",vertex.label.cex=1.3,layout=l*1,rescale=F)
legend("topleft", c("AMR","CMR", "STA","Clin"), pch=20,
       col=COLOR[1:4],  pt.cex=5,cex=1.5, bty="n", ncol=1)
legend("bottomleft", c("Neg-association","Pos-associatiom"), lty=1,lwd=5,
       col=c("darkgoldenrod3","gray56"), cex=1.5, bty="n", ncol=1)
dev.off()


colnames(nodes)[1]<-"Genes"
merge(edges,nodes,by="Genes")
write.csv(merge(edges,nodes,by="Genes"),file="Results/RNAseq/genes_associated_clin.csv")


#################################################################
## Study how non-coding genes are associated with coding-genes ##
################################################################

##separate between coding and non-coding
results_coding<-results_multi[which(results_multi$type_gene=="protein_coding"),]
results_noncoding<-results_multi[which(results_multi$type_gene!="protein_coding"),]

##norm data with coding
id_sign_coding<-match(results_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[id_sign_coding,]
rownames(norm_data_rlog_coding)<-results_coding$name

##norm data with non-coding
id_sign_noncoding<-match(results_noncoding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_noncoding<-t(norm_data_rlog[id_sign_noncoding,])
colnames(norm_data_rlog_noncoding)<-results_noncoding$name

###Find which non-coding are associated with coding genes
####Linear regression

##1. Find all the p-values
pvalue<-matrix(NA,ncol(norm_data_rlog_noncoding),nrow(norm_data_rlog_coding))
for (i in 1:ncol(norm_data_rlog_noncoding)){
  for (j in 1:nrow(norm_data_rlog_coding)){
    pvalue[i,j]<-coef(summary(lm(norm_data_rlog_noncoding[,i]~norm_data_rlog_coding[j,])))[2,4]
  }
}
pvalue_adj<-p.adjust(pvalue,method = "fdr")  
maxp<-max(pvalue_adj[which(pvalue_adj<0.05)])
grep(maxp, pvalue_adj) #position 1096
pvalue_array<-as.numeric(pvalue)
fdr_value<-pvalue_array[1096] #pvalue=0.01070372 ##This is the FDR correction < 0.05

##2. Build the matrix
lm_pvalue<-rep(NA,ncol(norm_data_rlog_coding))
lm_coef<-rep(NA,ncol(norm_data_rlog_coding))
lm_noncoding_p<-list()
lm_noncoding_c<-list()
for (i in 1:ncol(norm_data_rlog_noncoding)){
  for (j in 1:nrow(norm_data_rlog_coding)){
    lm_pvalue[j]<-coef(summary(lm(norm_data_rlog_noncoding[,i]~norm_data_rlog_coding[j,])))[2,4]
    lm_coef[j]<-coef(summary(lm(norm_data_rlog_noncoding[,i]~norm_data_rlog_coding[j,])))[2,1]
  }
  sign_pvalue<-lm_pvalue[which(lm_pvalue<fdr_value)]
  sign_coef<-lm_coef[which(lm_pvalue<fdr_value)]
  names(sign_pvalue)<-rownames(norm_data_rlog_coding)[which(lm_pvalue<fdr_value)]
  names(sign_coef)<-rownames(norm_data_rlog_coding)[which(lm_pvalue<fdr_value)]
  lm_noncoding_p[[i]]<-sign_pvalue
  lm_noncoding_c[[i]]<-sign_coef
}
names(lm_noncoding_p)<-colnames(norm_data_rlog_noncoding)
names(lm_noncoding_c)<-colnames(norm_data_rlog_noncoding)

noncoding_df<- data.frame(
  Genes_noncod = rep(names(lm_noncoding_p), lapply(lm_noncoding_p, length)), Genes_cod = unlist(lapply(lm_noncoding_p, names)),
  pvalue = unlist(lm_noncoding_p), coef = unlist(lm_noncoding_p))

join_non_coding<-inner_join(noncoding_df,results_multi,by=c("Genes_noncod"="name"))
join_non_coding<-join_non_coding[,c(1:4,10,14)]
colnames(join_non_coding)[5:6]<-c("type_noncod","cluster_noncod")
join_coding<-inner_join(join_non_coding,results_multi,by=c("Genes_cod"="name"))
join_coding<-join_coding[,c(1:6,16)]
colnames(join_coding)[7]<-c("cluster_cod")
noncoding_df<-join_coding
write.csv(noncoding_df,"Results/RNAseq/results_noncoding.csv")

###After looking for the GO biological terms in EnrichR
noncoding_df_GO<-read.csv("Results/RNAseq/results_noncoding_GO.csv")
table(noncoding_df_GO$GO.Biological.Process)

GO_terms <- group_by(noncoding_df_GO,GO.Biological.Process)
summarize(GO_terms)

##Build a co-expression matrix
cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}

###Find the significant genes in the normalize data
coding_genes<-unique(as.character(noncoding_df$Genes_cod))
coding_genes[26]<-"SEPT2"
non_coding_genes<-unique(as.character(noncoding_df$Genes_noncod))
genes_results<-c(coding_genes,non_coding_genes)
id<-match(genes_results,annotation$name)
norm_data_rlog_results<-norm_data_rlog[match(annotation$Ensemble_id[id],rownames(norm_data_rlog)),]
id_genes<-match(rownames(norm_data_rlog_results),annotation$Ensemble_id)
rownames(norm_data_rlog_results)<-annotation$name[id_genes]

###Build de similarity matrix
sim_matrix <- cordist(norm_data_rlog_results)
id_coding<-match(coding_genes,rownames(sim_matrix))
id_noncoding<-match(non_coding_genes,rownames(sim_matrix))
sim_matrix_2<-sim_matrix[id_coding,id_noncoding]
  
###Find which cluster is up-regulated
results_multi$name<-as.character(results_multi$name)
results_multi$name[19]<-c("SEPT2")
cluster_coding<-results_multi$cluster[match(rownames(sim_matrix_2),results_multi$name)]
cluster_noncoding<-results_multi$cluster[match(colnames(sim_matrix_2),results_multi$name)]
COLOR = brewer.pal(4,"Pastel1")
ann_colors = list(cluster_coding = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]),
                  cluster_noncoding  = c("AMR" = COLOR [1] ,"CMR" = COLOR[2], "STA" = COLOR[3]))
names(cluster_coding)<-rownames(sim_matrix_2)
names(cluster_noncoding)<-colnames(sim_matrix_2)

tiff("Results/RNAseq/SimilarityMatrix.tiff",res=300,h=3000,w=4000)
pheatmap(sim_matrix_2,border_color=F,annotation_colors = ann_colors,
         annotation_row = as.data.frame(cluster_coding),annotation_col=as.data.frame(cluster_noncoding))
dev.off()


################################
## Match with microarray data ##
###############################
geneMicroarray<-read.csv("Data/MicroaarayGeneList.csv")
match(results_multi$name,geneMicroarray$GeneSymbol)
match(results_AMR_CMR$name,geneMicroarray$GeneSymbol)
match(results_AMR_STA$name,geneMicroarray$GeneSymbol)
geneMicroarray$GeneSymbol[na.omit(match(results_REJ_STA$name,geneMicroarray$GeneSymbol))]

geneKsort<-read.csv("Data/kSORT_list.csv")
match(results_REJ_STA$name,geneKsort$Genes)
match(results_STA_REJ_DEseq$name,geneKsort$Genes)

geneKurian<-read.csv("Data/kurian_genelist.csv")
geneKurian$Symbol[na.omit(match(results_REJ_STA$name,geneKurian$Symbol))]

#############################
## Match with ExomeSeq #####
############################
geneExomeSeq<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",header=T,sep="\t")
match(results_multi$name,geneExomeSeq$Gene.refGene)
match(results_AMR_CMR$name,geneExomeSeq$Gene.refGene)
match(results_AMR_STA$name,geneExomeSeq$Gene.refGene)
##None match with this list