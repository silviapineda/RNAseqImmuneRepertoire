###Read the file in gtf 
fh = open("gencode.v27lift37.annotation.gtf")
fout = open("gencode.v27.annotation_genes.txt",'w')
count = 0
for line in fh:
	line=line.rstrip()
	if line.startswith("##"):
		fout.write( line + '\n')
		continue
  	words = line.split()
  	#print words
	if ("gene") in words[2]:
	    fout.write( words[0] + ';' + words[3] + ';' + words[4] + ';' + words[9]+words[11]+words[13] + '\n')
	    count=count+1
print count
print "Done"

