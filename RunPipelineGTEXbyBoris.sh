echo
echo SRR$srrID
#/opt/sratoolkit.2.5.2/bin/prefetch --max-size 1000000000 -a "/home/oskotsky/.aspera/connect/bin/ascp|/home/oskotsky/.aspera/connect/etc/asperaweb_id_dsa.openssh" SRR$srrID
/opt/sratoolkit.2.5.2/bin/fastq-dump --split-files $sraDir/SRR$srrID.sra -O $fastqDir
ext1="_1.fastq"
ext2="_2.fastq"
fastqFilePath1=SRR$srrID$ext1
fastqFilePath2=SRR$srrID$ext2


###### Pipeline to run MIXCR
export PATH=$PATH:/home/pinedasans/programs/mixcr-2.1.3/
##Align the two fastq files
mixcr align -p rna-seq -s hsa -f -r {SRR$srrID}_aligments_report.txt -t 20 -OallowPartialAlignments=true SRR$srrID$ext1 SRR$srrID$ext2 {SRR$srrID}_alignments.vdjca


done
rm $sraDir/SRR$srrID.sra
rm $fastqDir/$fastqFilePath1
rm $fastqDir/$fastqFilePath2
echo SRR$srrID Done...
