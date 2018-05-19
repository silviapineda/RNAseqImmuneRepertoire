
###Read the file in gtf 
fh = open("gencode.v26.annotation.gtf")
fout = open("gencode.v26.annotation_coding_genes.txt",'w')
count = 0
for line in fh:
	line=line.rstrip()
	if line.startswith("##"):
		fout.write( line + '\n')
		continue
  	words = line.split()
	if ("gene") not in words[2]: continue
	if ("protein_coding") in line:
		fout.write( words[0] + ';' + words[3] + ';' + words[4] + ';' + words[9]+words[11]+words[17] + '\n')
		count=count+1
print count
print "Done"


