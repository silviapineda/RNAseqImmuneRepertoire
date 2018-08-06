#Work Directory
setwd("/Users/Pinedasans/Catalyst/Data/ExomeSeq/")

###Reading the Total genotypes data frame
load("ExomeSeqVCF_SNPs.Rdata")
demographics<-read.table("/Users/Pinedasans/Catalyst/Data/Demographics.txt",sep="\t",header=T)
non.list<- seq(1,56,2)
clin<-demographics$phenotype[non.list]

###SNPs in the AK9 gene
variants_AK9<-df_joint_qc[grep("AK9",df_joint_qc$Gene.refGene),c(17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71)]

p<-NULL
for(i in 60:nrow(variants_AK9)){
  print(i)
  p[i]<-fisher.test(table(variants_AK9[i,],clin))$p.value
}
