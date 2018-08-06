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

working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

###Load the data
load("Data/norm_data_rlog_filter.Rdata")

results_multi<-read.csv("Results/RNAseq/genes.enet.multinomial.csv")
results_AMR_STA<-read.csv("Results/RNAseq/results_STA_AMR_ENET.csv")
results_AMR_CMR<-read.csv("Results/RNAseq/results_CMR_AMR_ENET.csv")

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
id_genes<-match(results$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_results<-norm_data_rlog[id_genes,]

#put the names of the genes
id_genes<-match(rownames(norm_data_rlog_results),results$Ensemble_id)
rownames(norm_data_rlog_results)<-results$name[id_genes]

clin_vars<-clinData[,c(5:17)]
colnames(clin_vars)<-c("ag","ai","at","ti","ptc","av","aah","cg","ci","ct","cv","cm","C4d")
rownames(clin_vars)<-clinData$Individual_id
####Correlation analysis
cor<-rep(NA,ncol(clin_vars))
p<-rep(NA,ncol(clin_vars))
cor_network<-list()
for (i in 1:nrow(norm_data_rlog_results)){
  for (j in 1:ncol(clin_vars)){
    p[j]<-cor.test(norm_data_rlog_results[i,],clin_vars[,j],use="na.or.complete")$p.value
    cor[j]<-cor(norm_data_rlog_results[i,],clin_vars[,j],use="na.or.complete","spearman")
  }
  sign<-cor[which(p<0.01)]
  names(sign)<-colnames(clin_vars)[which(p<0.01)]
  cor_network[[i]]<-sign
}
names(cor_network)<-rownames(norm_data_rlog_results)

###Multivarate analysis per gene
p_network<-list()
for(i in 1:nrow(norm_data_rlog_results)){
  lm_results<-coef(summary(lm(norm_data_rlog_results[i,]~ag+ai+at+ti+ptc+av+aah+cg+ci+ct+cv+cm+C4d,
                              data = clin_vars)))[,4]
  sign<-lm_results[which(lm_results<0.01)]
  p_network[[i]]<-sign[-1]
    
}

names(p_network)<-rownames(norm_data_rlog_results)

####Linear regression
lm_pvalue<-rep(NA,ncol(clin_vars))
lm_coef<-rep(NA,ncol(clin_vars))
lm_network_p<-list()
lm_network_c<-list()
for (i in 1:nrow(norm_data_rlog_results)){
  for (j in 1:ncol(clin_vars)){
    lm_pvalue[j]<-coef(summary(lm(norm_data_rlog_results[i,]~clin_vars[,j], data = clin_vars)))[2,4]
    lm_coef[j]<-coef(summary(lm(norm_data_rlog_results[i,]~clin_vars[,j], data = clin_vars)))[2,1]
  }
  sign_pvalue<-lm_pvalue[which(lm_pvalue<0.01)]
  sign_coef<-lm_coef[which(lm_pvalue<0.01)]
  names(sign_pvalue)<-colnames(clin_vars)[which(lm_pvalue<0.01)]
  names(sign_coef)<-colnames(clin_vars)[which(lm_pvalue<0.01)]
  lm_network_p[[i]]<-sign_pvalue
  lm_network_c[[i]]<-sign_coef
}
names(lm_network_p)<-rownames(norm_data_rlog_results)
names(lm_network_c)<-rownames(norm_data_rlog_results)

##Build the network
# ##For p-values
# edges<- data.frame(
#   Genes = rep(names(p_network), lapply(p_network, length)),Clin = unlist(lapply(p_network, names)),
#   value = abs(log(unlist(p_network),10)))
# ##For correlation values
# edges<- data.frame(
#   Genes = rep(names(cor_network), lapply(cor_network, length)),Clin = unlist(lapply(cor_network,names)),
#   value = abs(unlist(cor_network)))

##For lm
edges<- data.frame(
  Genes = rep(names(lm_network_p), lapply(lm_network_p, length)), Clin = unlist(lapply(lm_network_p, names)),
  pvalue = abs(log(unlist(lm_network_p),10)), coef = unlist(lm_network_c))

edges$color<-ifelse(edges$coef<0,"darkgoldenrod3","gray56")

#Unique Genes with edges
Genes<-as.character(unique(edges$Genes))
#All genes significant in the DE analysis
#Genes<-as.character(results$name)
id_AMR<-na.omit(match(results$name[which(results$cluster=="AMR")],Genes))
id_CMR<-na.omit(match(results$name[which(results$cluster=="CMR")],Genes))
id_STA<-na.omit(match(results$name[which(results$cluster=="STA")],Genes))

nodes<-c(Genes[id_AMR],Genes[id_CMR],Genes[id_STA],colnames(clin_vars))
COLOR = brewer.pal(4,"Pastel1")
nodes<-cbind(nodes,c(rep(COLOR[1],length(id_AMR)),rep(COLOR[2],length(id_CMR)),
                     rep(COLOR[3],length(id_STA)),rep(COLOR[4],ncol(clin_vars))))
nodes<-data.frame(cbind(nodes,c(rep(5,length(Genes)),rep(10,ncol(clin_vars)))))
colnames(nodes)<-c("names","color","size")
nodes$size<-as.numeric(as.character(nodes$size))



#Make the graph and plot
net <- graph_from_data_frame(d = edges, vertices = nodes,directed = F)

E(net)$width<-E(net)$pvalue
tiff(paste("Results/RNAseq/network_lm_01.tiff",sep=""),res=300,h=5000,w=5000)
l <- layout_with_fr(net)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net,vertex.label.color="black",vertex.label.cex=1.3,layout=l*1)
legend("topleft", c("AMR","CMR", "STA","Clin"), pch=20,
       col=COLOR[1:4],  pt.cex=5,cex=1.5, bty="n", ncol=1)
legend("bottomleft", c("Neg-association","Pos-associatiom"), lty=1,lwd=5,
       col=c("darkgoldenrod3","gray56"), cex=1.5, bty="n", ncol=1)
dev.off()

##Try other ways of layout
e <- get.edgelist(net,names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net))
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net),
            area=8*(vcount(net)^2),repulse.rad=(vcount(net)^3.1))
plot(net,layout=l,vertex.size=4,vertex.label=NA)


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
lm_pvalue<-rep(NA,ncol(norm_data_rlog_coding))
lm_coef<-rep(NA,ncol(norm_data_rlog_coding))
lm_noncoding_p<-list()
lm_noncoding_c<-list()
for (i in 1:ncol(norm_data_rlog_noncoding)){
  for (j in 1:nrow(norm_data_rlog_coding)){
    lm_pvalue[j]<-coef(summary(lm(norm_data_rlog_noncoding[,i]~norm_data_rlog_coding[j,])))[2,4]
    lm_coef[j]<-coef(summary(lm(norm_data_rlog_noncoding[,i]~norm_data_rlog_coding[j,])))[2,1]
  }
  sign_pvalue<-lm_pvalue[which(lm_pvalue<0.01)]
  sign_coef<-lm_coef[which(lm_pvalue<0.01)]
  names(sign_pvalue)<-rownames(norm_data_rlog_coding)[which(lm_pvalue<0.01)]
  names(sign_coef)<-rownames(norm_data_rlog_coding)[which(lm_pvalue<0.01)]
  lm_noncoding_p[[i]]<-sign_pvalue
  lm_noncoding_c[[i]]<-sign_coef
}
names(lm_noncoding_p)<-colnames(norm_data_rlog_noncoding)
names(lm_noncoding_c)<-colnames(norm_data_rlog_noncoding)


##For lm
edges<- data.frame(
  Genes = rep(names(lm_noncoding_p), lapply(lm_noncoding_p, length)), coding = unlist(lapply(lm_noncoding_p, names)),
  pvalue = abs(log(unlist(lm_noncoding_p),10)), coef = unlist(lm_noncoding_c))

edges$color<-ifelse(edges$coef<0,"darkgoldenrod3","gray56")

#Unique Genes with edges
Genes<-c(as.character(unique(edges$Genes)),as.character(edges$coding))
#All genes significant in the DE analysis
#Genes<-as.character(results$name)
id_AMR<-na.omit(match(results_multi$name[which(results_multi$cluster=="AMR")],Genes))
id_CMR<-na.omit(match(results_multi$name[which(results_multi$cluster=="CMR")],Genes))
id_STA<-na.omit(match(results_multi$name[which(results_multi$cluster=="STA")],Genes))

nodes<-c(Genes[id_AMR],Genes[id_CMR],Genes[id_STA],colnames(norm_data_rlog_noncoding))
COLOR = brewer.pal(4,"Pastel1")
nodes<-cbind(nodes,c(rep(COLOR[1],length(id_AMR)),rep(COLOR[2],length(id_CMR)),
                     rep(COLOR[3],length(id_STA)),rep(COLOR[4],ncol(norm_data_rlog_noncoding))))
nodes<-data.frame(cbind(nodes,c(rep(5,length(Genes)),rep(5,ncol(norm_data_rlog_noncoding)))))
colnames(nodes)<-c("names","color","size")
nodes$size<-as.numeric(as.character(nodes$size))

#Make the graph and plot
net <- graph_from_data_frame(d = edges, vertices = nodes,directed = F)

E(net)$width<-E(net)$pvalue
tiff(paste("Results/RNAseq/network_coding-noncoding.tiff",sep=""),res=300,h=5000,w=5000)
l <- layout_with_fr(net)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net,vertex.label.color="black",vertex.label.cex=1.3,layout=l*1)
legend("topleft", c("AMR","CMR", "STA","Clin"), pch=20,
       col=COLOR[1:4],  pt.cex=5,cex=1.5, bty="n", ncol=1)
legend("bottomleft", c("Neg-association","Pos-associatiom"), lty=1,lwd=5,
       col=c("darkgoldenrod3","gray56"), cex=1.5, bty="n", ncol=1)
dev.off()



