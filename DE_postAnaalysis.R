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


working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)

###Load the data
load("Data/norm_data_rlog_filter.Rdata")

results<-read.csv("Results/RNAseq/genes.enet.multinomial.csv")

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
id_filter_genes<-match(rownames(norm_data_rlog),annotation$Ensemble_id)
annotation_filter<-annotation[id_filter_genes,] #26,545

###Non-coding
results$geneType<-ifelse(results$type_gene=="protein_coding","coding","non-coding")
tab<-table(results$geneType,results$cluster)
chisq.test(tab) #p-value = 0.225

###lincRNA
results_lncRNA<-results[which(results$type_gene=="protein_coding" | results$type_gene=="lincRNA"),]
tab<-table(results_lncRNA$geneType,results_lncRNA$cluster)
fisher.test(tab) #p-value = 0.5


tab_total<-table(annotation_filter$type_gene)
dim(annotation_filter)

#AMR
counts=matrix(data=c(32,18,(15420-32),(11125-18)),nrow=2)
chisq.test(counts) #p-value =0.5
#CMR
counts=matrix(data=c(16,6,(15420-16),(11125-6)),nrow=2)
chisq.test(counts) #p-value =0.2
#STA
counts=matrix(data=c(15,15,(15420-15),(11125-15)),nrow=2)
chisq.test(counts) #p-value =0.5

###Correlation with the clinical variables
id_genes<-match(results$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_results<-norm_data_rlog[id_genes,]

clin_vars<-clinData[,c(5:17)]

cor<-matrix(NA,nrow(norm_data_rlog_results),ncol(clin_vars))
for (i in 1:nrow(norm_data_rlog_results)){
  for (j in 1:ncol(clin_vars)){
    cor[i,j]<-cor(norm_data_rlog_results[i,],clin_vars[,j],use="na.or.complete")
  }
}

colnames(cor)<-colnames(clin_vars)
rownames(cor)<-rownames(norm_data_rlog_results)
