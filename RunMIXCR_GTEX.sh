#!/bin/bash
##Run the MIXCR program to find immuneSEq from RNAseq in GTEX
##https://mixcr.readthedocs.io/en/latest/rnaseq.html#ref-rna-seq

###### Pipeline to run MIXCR
export PATH=$PATH:/home/pinedasans/programs/mixcr-2.1.5/

## 1. Align sequencing reads
##Run by Boris in the server
mixcr align -p rna-seq -s hsa -f -r {SRR$srrID}_aligments_report.txt -t 20 -OallowPartialAlignments=true SRR$srrID$ext1 SRR$srrID$ext2 {SRR$srrID}_alignments.vdjca


#2. Partial alignments
for i in $(ls *alignments.vdjca | rev | cut -c 17- | rev | uniq)

do
echo $i 

mixcr assemblePartial -f -r ${i}aligments_report_rescued_1.txt ${i}alignments.vdjca ${i}alignments_rescued_1.vdjca
mixcr assemblePartial -f -r ${i}aligments_report_rescued_2.txt ${i}alignments_rescued_1.vdjca ${i}alignments_rescued_2.vdjca

#mixcr assemblePartial -f -r {SRR1305191}_aligments_report_rescued_1.txt {SRR1305191}_alignments.vdjca {SRR1305191}_alignments_rescued_1.vdjca
#mixcr assemblePartial -f -r {SRR1305191}_aligments_report_rescued_2.txt {SRR1305191}_alignments_rescued_1.vdjca {SRR1305191}_alignments_rescued_2.vdjca

done

##3. Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr extendAlignments -f -r ${i}aligments_report_extended.txt  ${i}alignments_rescued_2.vdjca ${i}alignments_extended.vdjca

#mixcr extendAlignments -f -r {SRR1305191}_aligments_report_extended.txt  {SRR1305191}_alignments_rescued_2.vdjca {SRR1305191}_alignments_extended.vdjca

done

##4. Assemble the clonotypes
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
mixcr assemble -r ${i}clonotypes_report.txt -i ${i}index_file ${i}alignments_extended.vdjca ${i}output.clns

#mixcr assemble -f -r {SRR1305191}_clonotypes_report.txt -i {SRR1305191}_index_file {SRR1305191}_alignments_rescued_2.vdjca {SRR1305191}_output.clns

done

##5. Extract Alignments
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
#mixcr exportAlignments -readId -sequence -cloneId ${i}index_file -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene ${i}alignments_extended.vdjca ${i}alignments.txt
mixcr exportAlignments -readId -sequence -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene ${i}alignments_extended.vdjca ${i}alignments.txt
#mixcr exportAlignments -f -readId -sequence -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene {SRR1305191}_alignments_extended.vdjca {SRR1305191}_alignments.txt

done

##6. Extract Clones
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)
do
echo $i
mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${i}output.clns ${i}clones.txt

#mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3  {SRR3478950}_output.clns {SRR3478950}_clones.txt
done