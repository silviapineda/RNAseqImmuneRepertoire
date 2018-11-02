working_directory<-"/Users/Pinedasans/ImmuneRep_RNAseq/"
setwd(working_directory)


# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000717","nucleotide-excision repair, DNA duplex unwinding",0.001,4.5768,0.645,0.000,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:1901292","nucleoside phosphate catabolic process",0.187,2.2265,0.698,0.200,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0006165","nucleoside diphosphate phosphorylation",0.600,2.1177,0.712,0.532,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:1903506","regulation of nucleic acid-templated transcription",9.965,2.8333,0.589,0.537,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:1904153","negative regulation of retrograde protein transport, ER to cytosol",0.002,2.1440,0.568,0.676,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0006294","nucleotide-excision repair, preincision complex assembly",0.001,4.3420,0.678,0.494,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:1900101","regulation of endoplasmic reticulum unfolded protein response",0.007,2.3019,0.693,0.326,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0009223","pyrimidine deoxyribonucleotide catabolic process",0.001,1.9586,0.745,0.638,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0031053","primary miRNA processing",0.003,2.0571,0.606,0.500,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0006359","regulation of transcription from RNA polymerase III promoter",0.033,1.9416,0.700,0.270,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0032508","DNA duplex unwinding",0.590,4.0467,0.775,0.618,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0006396","RNA processing",3.210,2.7203,0.711,0.164,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0016071","mRNA metabolic process",0.798,2.5643,0.734,0.373,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0070911","global genome nucleotide-excision repair",0.003,4.2581,0.718,0.529,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0031106","septin ring organization",0.009,2.0474,0.790,0.311,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0006397","mRNA processing",0.561,2.3951,0.681,0.358,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0046033","AMP metabolic process",0.124,2.1687,0.726,0.468,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0044806","G-quadruplex DNA unwinding",0.003,2.5533,0.808,0.649,"nucleotide-excision repair, DNA duplex unwinding"),
                     c("GO:0044788","modulation by host of viral process",0.004,2.0474,0.780,0.000,"modulation by host of viral process"),
                     c("GO:1902188","positive regulation of viral release from host cell",0.003,1.9230,0.775,0.624,"modulation by host of viral process"),
                     c("GO:0055076","transition metal ion homeostasis",0.197,1.9055,0.817,0.402,"modulation by host of viral process"),
                     c("GO:0050657","nucleic acid transport",0.100,2.1950,0.891,0.000,"nucleic acid transport"),
                     c("GO:0070972","protein localization to endoplasmic reticulum",0.187,1.9527,0.873,0.209,"nucleic acid transport"),
                     c("GO:1902224","ketone body metabolic process",0.013,2.2606,0.910,0.022,"ketone body metabolism"),
                     c("GO:0046951","ketone body biosynthetic process",0.002,2.3476,0.832,0.027,"ketone body biosynthesis"),
                     c("GO:1901570","fatty acid derivative biosynthetic process",0.009,1.9055,0.857,0.195,"ketone body biosynthesis"),
                     c("GO:0060972","left/right pattern formation",0.006,2.1687,0.915,0.039,"left/right pattern formation"),
                     c("GO:0019184","nonribosomal peptide biosynthetic process",0.076,2.1883,0.799,0.093,"nonribosomal peptide biosynthesis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
tiff( file="Results/RNAseq/revigo_treemap_STA.tiff", res=300,width=2600, height=1900 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Gene Ontology treemap for non-coding genes up-regulated in STA",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
