#!/bin/bash
##Run the MIXCR program to find immuneSEq from RNAseq
export PATH=$PATH:/home/pinedasans/programs/mixcr-2.1.5/

#First alignement
for i in $(ls ~/ExomeSeq/FASTQ_secondrun/*.fastq.gz | rev | cut -c 17- | rev | uniq)
do
echo $i
mixcr align -OvParameters.geneFeatureToAlign=VGeneWithP -s hsa -f -r ${i}_aligments_report.txt -t 20 -OallowPartialAlignments=true ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz ${i}_alignments.vdjca

echo ${i}.fastq
done

#Partial alignments
for i in $(ls *alignments.vdjca | rev | cut -c 18- | rev | uniq)

do
echo $i

mixcr assemblePartial -f -r ${i}_alignments_report_rescued_1.txt ${i}_alignments.vdjca ${i}_alignments_rescued_1.vdjca
mixcr assemblePartial -f -r ${i}_alignments_report_rescued_2.txt ${i}_alignments_rescued_1.vdjca ${i}_alignments_rescued_2.vdjca

done

##Extend alignments
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 28- | rev | uniq)

do
echo $i
mixcr extendAlignments -f -r ${i}_alignments_report_extended.txt  ${i}_alignments_rescued_2.vdjca ${i}_alignments_extended.vdjca

done

##Assemble the clonotypes
for i in $(ls *alignments_extended.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr assemble -f -r ${i}_clonotypes_report.txt  -ObadQualityThreshold=15 ${i}_alignments_extended.vdjca ${i}_output.clns

done

###Try to run the merge alignments


##Extract Alignments
for i in $(ls *alignments_extended.vdjca | rev | cut -c 27- | rev | uniq)

#Run with the options to export
do
echo $i
mixcr exportAlignments -f -readId -sequence -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -defaultAnchorPoints -nFeature CDR3 -aaFeature CDR3 ${i}_al
ignments_extended.vdjca ${i}_alignments.txt

done

##Extract Clones
##6. Extract Clones
for i in $(ls *alignments_extended.vdjca | rev | cut -c 27- | rev | uniq)
do
echo $i
mixcr exportClones -f -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${i}_output.clns ${i}_clones.txt

done