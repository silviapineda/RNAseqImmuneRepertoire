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
revigo.data <- rbind(c("GO:1903311","regulation of mRNA metabolic process",0.044,2.6202,0.515,0.000,"regulation of mRNA metabolism"),
                     c("GO:0006165","nucleoside diphosphate phosphorylation",0.600,2.2018,0.659,0.516,"regulation of mRNA metabolism"),
                     c("GO:0045191","regulation of isotype switching",0.006,2.3196,0.343,0.538,"regulation of mRNA metabolism"),
                     c("GO:0006306","DNA methylation",0.190,1.8788,0.636,0.173,"regulation of mRNA metabolism"),
                     c("GO:0032785","negative regulation of DNA-templated transcription, elongation",0.005,2.1561,0.551,0.235,"regulation of mRNA metabolism"),
                     c("GO:0035413","positive regulation of catenin import into nucleus",0.002,2.2686,0.612,0.143,"regulation of mRNA metabolism"),
                     c("GO:2000095","regulation of Wnt signaling pathway, planar cell polarity pathway",0.003,2.3196,0.545,0.148,"regulation of mRNA metabolism"),
                     c("GO:0006355","regulation of transcription, DNA-templated",9.917,3.3537,0.484,0.413,"regulation of mRNA metabolism"),
                     c("GO:0046033","AMP metabolic process",0.124,2.2528,0.654,0.135,"regulation of mRNA metabolism"),
                     c("GO:0060972","left/right pattern formation",0.006,2.2528,0.760,0.402,"regulation of mRNA metabolism"),
                     c("GO:0050684","regulation of mRNA processing",0.035,2.3621,0.473,0.611,"regulation of mRNA metabolism"),
                     c("GO:0006915","apoptotic process",0.406,1.4634,0.798,0.153,"regulation of mRNA metabolism"),
                     c("GO:0010559","regulation of glycoprotein biosynthetic process",0.008,2.1561,0.622,0.320,"regulation of mRNA metabolism"),
                     c("GO:1903513","endoplasmic reticulum to cytosol transport",0.011,1.9266,0.857,0.412,"regulation of mRNA metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
tiff( file="Results/RNAseq/revigo_treemap_AMR.tiff", res=300,width=2600, height=1900 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Gene Ontology treemap for non-coding genes up-regulated in AMR",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
