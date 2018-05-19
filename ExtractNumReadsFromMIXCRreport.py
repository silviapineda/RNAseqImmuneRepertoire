import glob

###Read the file in gtf 
filenames = glob.glob("*_aligments_report.txt")
fout = open("/home/pinedasans/ImmuneRep_RNAseq/total_reads_GTEX.txt",'w')
for f in filenames:
    print(f)
    sample = f[1:11]
    print(sample)
    fh = open(f)
    for line in fh:
        line=line.rstrip()
        if line.startswith("Total sequencing reads"):
            words = line.split()
            print(words[3])
            fout.write(sample + ';' + words[3] + '\n')
            continue
print("Done")