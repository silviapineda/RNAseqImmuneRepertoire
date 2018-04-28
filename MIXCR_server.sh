#!/bin/bash
##Run the MIXCR program to find immuneSEq from RNAseq

#First alignement
for i in $(ls *.fastq.gz | rev | cut -c 12- | rev | uniq)

do
echo $i
mixcr align -p rna-seq -s hsa -f -r ${i}aligments_report.txt -t 20 -OallowPartialAlignments=true ${i}R1.fastq.gz ${i}R2.fastq.gz ${i}alignments.vdjca

echo ${i}R1.fastq.gz
echo ${i}R2.fastq.gz
done

#Partial alignments
for i in $(ls *alignments.vdjca | rev | cut -c 17- | rev | uniq)

do
echo $i 

mixcr assemblePartial -f -r ${i}aligments_report_rescued_1.txt ${i}alignments.vdjca ${i}alignments_rescued_1.vdjca
mixcr assemblePartial -f -r ${i}aligments_report_rescued_2.txt ${i}alignments_rescued_1.vdjca ${i}alignments_rescued_2.vdjca

done

##Extend alignments 
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr extendAlignments -f -r ${i}aligments_report_extended.txt  ${i}alignments_rescued_2.vdjca ${i}alignments_extended.vdjca

done

##Assemble the clonotypes
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
mixcr assemble -r ${i}clonotypes_report.txt -i ${i}index_file -ObadQualityThreshold=15 ${i}alignments_extended.vdjca ${i}output.clns

done

###Try to run the merge alignments


##Extract Alignments
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

##Run with the options to export
do
echo $i
mixcr exportAlignments -f -readId -sequence -cloneId ${i}index_file -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -defaultAnchorPoints -nFeature CDR3 -aaFeature CDR3 ${i}alignments_extended.vdjca ${i}alignments.txt

done

##Extract Clones
##6. Extract Clones
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)
do
echo $i
mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${i}output.clns ${i}clones.txt

done