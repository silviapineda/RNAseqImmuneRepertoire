
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
revigo.data <- rbind(c("GO:0046033","AMP metabolic process",0.124,2.9208,0.523,0.000,"AMP metabolism"),
                     c("GO:0006165","nucleoside diphosphate phosphorylation",0.600,2.8700,0.563,0.516,"AMP metabolism"),
                     c("GO:0010907","positive regulation of glucose metabolic process",0.009,2.2229,0.635,0.200,"AMP metabolism"),
                     c("GO:0006754","ATP biosynthetic process",0.432,2.3061,0.468,0.639,"AMP metabolism"),
                     c("GO:0006355","regulation of transcription, DNA-templated",9.917,4.0155,0.472,0.339,"AMP metabolism"),
                     c("GO:0010559","regulation of glycoprotein biosynthetic process",0.008,2.3019,0.533,0.320,"AMP metabolism"),
                     c("GO:0006768","biotin metabolic process",0.081,2.6202,0.656,0.232,"AMP metabolism"),
                     c("GO:0000338","protein deneddylation",0.017,2.6580,0.668,0.039,"protein deneddylation"),
                     c("GO:0006306","DNA methylation",0.190,2.4208,0.461,0.280,"protein deneddylation"),
                     c("GO:0006296","nucleotide-excision repair, DNA incision, 5'-to lesion",0.003,1.7249,0.626,0.324,"protein deneddylation"),
                     c("GO:0043968","histone H2A acetylation",0.010,2.4690,0.673,0.264,"protein deneddylation"),
                     c("GO:0042990","regulation of transcription factor import into nucleus",0.020,2.3986,0.747,0.054,"regulation of transcription factor import into nucleus"),
                     c("GO:0002029","desensitization of G-protein coupled receptor protein signaling pathway",0.002,2.3476,0.724,0.139,"regulation of transcription factor import into nucleus"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
tiff( file="Results/RNAseq/revigo_treemap_CMR.tiff", res=300,width=2600, height=1900 ) # width and height are in inches
# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Gene Ontology treemap for non-coding genes up-regulated in TCMR",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
