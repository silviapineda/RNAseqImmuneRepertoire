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

done

##3. Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr extendAlignments -f -r ${i}aligments_report_extended.txt  ${i}alignments_rescued_2.vdjca ${i}alignments_extended.vdjca

done

##4. Assemble the clonotypes
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
mixcr assemble -r ${i}clonotypes_report.txt -i ${i}index_file ${i}alignments_extended.vdjca ${i}output.clns

done

##5. Extract Alignments
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
mixcr exportAlignments -readId -sequence -cloneId ${i}index_file -vGenes -dGenes -jGenes ${i}alignments_extended.vdjca ${i}alignments.txt

done

##6. Extract Clones
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)
do
echo $i
mixcr exportClones -readId -sequence -cloneId ${i}index_file ${i}output.clns ${i}clones.txt
done