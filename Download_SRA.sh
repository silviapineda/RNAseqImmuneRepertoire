for F in $(cat /home/pinedasans/ImmuneRep_RNAseq/VAL/SRR_Acc_List.txt) ; do
  echo $F
    echo SRR$srrID
    ./sratoolkit.2.9.0-ubuntu64/bin/fastq-dump $F

done
