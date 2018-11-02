

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
revigo.data <- rbind(c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.095,1.3510,0.971,0.000,"homophilic cell adhesion via plasma membrane adhesion molecules"),
                     c("GO:0007498","mesoderm development",0.031,1.8154,0.924,0.000,"mesoderm development"),
                     c("GO:0070228","regulation of lymphocyte apoptotic process",0.010,2.3107,0.806,0.000,"regulation of lymphocyte apoptotic process"),
                     c("GO:0000281","mitotic cytokinesis",0.070,1.3851,0.872,0.114,"regulation of lymphocyte apoptotic process"),
                     c("GO:0002718","regulation of cytokine production involved in immune response",0.013,2.3107,0.626,0.146,"regulation of lymphocyte apoptotic process"),
                     c("GO:0040036","regulation of fibroblast growth factor receptor signaling pathway",0.007,1.8355,0.589,0.582,"regulation of lymphocyte apoptotic process"),
                     c("GO:0000018","regulation of DNA recombination",0.045,1.6281,0.728,0.674,"regulation of lymphocyte apoptotic process"),
                     c("GO:0019932","second-messenger-mediated signaling",0.079,1.3643,0.682,0.337,"regulation of lymphocyte apoptotic process"),
                     c("GO:0048169","regulation of long-term neuronal synaptic plasticity",0.004,2.0773,0.770,0.140,"regulation of lymphocyte apoptotic process"),
                     c("GO:0097306","cellular response to alcohol",0.021,1.6965,0.698,0.564,"regulation of lymphocyte apoptotic process"),
                     c("GO:0071241","cellular response to inorganic substance",0.035,1.9020,0.707,0.600,"regulation of lymphocyte apoptotic process"),
                     c("GO:0070972","protein localization to endoplasmic reticulum",0.187,2.0106,0.884,0.559,"regulation of lymphocyte apoptotic process"),
                     c("GO:0007224","smoothened signaling pathway",0.038,1.6678,0.689,0.278,"regulation of lymphocyte apoptotic process"),
                     c("GO:0032418","lysosome localization",0.010,1.6819,0.907,0.176,"regulation of lymphocyte apoptotic process"),
                     c("GO:0030520","intracellular estrogen receptor signaling pathway",0.010,1.8566,0.633,0.458,"regulation of lymphocyte apoptotic process"),
                     c("GO:0070201","regulation of establishment of protein localization",0.165,2.2528,0.752,0.172,"regulation of lymphocyte apoptotic process"),
                     c("GO:0045579","positive regulation of B cell differentiation",0.003,2.2528,0.699,0.614,"regulation of lymphocyte apoptotic process"),
                     c("GO:0006953","acute-phase response",0.006,1.8787,0.774,0.197,"regulation of lymphocyte apoptotic process"),
                     c("GO:0000725","recombinational repair",0.152,1.3076,0.663,0.362,"regulation of lymphocyte apoptotic process"),
                     c("GO:0071478","cellular response to radiation",0.062,1.5274,0.739,0.229,"regulation of lymphocyte apoptotic process"),
                     c("GO:0072711","cellular response to hydroxyurea",0.001,2.2528,0.717,0.184,"regulation of lymphocyte apoptotic process"),
                     c("GO:0072710","response to hydroxyurea",0.002,2.2018,0.744,0.571,"regulation of lymphocyte apoptotic process"),
                     c("GO:0019722","calcium-mediated signaling",0.040,1.3135,0.692,0.413,"regulation of lymphocyte apoptotic process"),
                     c("GO:2000781","positive regulation of double-strand break repair",0.004,1.9529,0.553,0.466,"regulation of lymphocyte apoptotic process"),
                     c("GO:0031106","septin ring organization",0.009,1.9020,0.879,0.018,"septin ring organization"),
                     c("GO:0031346","positive regulation of cell projection organization",0.056,1.3576,0.700,0.290,"septin ring organization"),
                     c("GO:0043254","regulation of protein complex assembly",0.198,1.2959,0.741,0.595,"septin ring organization"),
                     c("GO:0016579","protein deubiquitination",0.195,1.8659,0.816,0.021,"protein deubiquitination"),
                     c("GO:0038083","peptidyl-tyrosine autophosphorylation",0.011,1.4985,0.813,0.478,"protein deubiquitination"),
                     c("GO:0018212","peptidyl-tyrosine modification",0.215,1.4461,0.791,0.577,"protein deubiquitination"),
                     c("GO:0043966","histone H3 acetylation",0.027,1.5175,0.760,0.327,"protein deubiquitination"),
                     c("GO:0009313","oligosaccharide catabolic process",0.019,2.1562,0.801,0.046,"oligosaccharide catabolism"),
                     c("GO:0009100","glycoprotein metabolic process",0.356,1.4715,0.863,0.466,"oligosaccharide catabolism"),
                     c("GO:1901136","carbohydrate derivative catabolic process",0.423,2.0773,0.817,0.418,"oligosaccharide catabolism"),
                     c("GO:0016052","carbohydrate catabolic process",1.078,1.6035,0.802,0.587,"oligosaccharide catabolism"),
                     c("GO:0016998","cell wall macromolecule catabolic process",0.059,2.0427,0.799,0.448,"oligosaccharide catabolism"),
                     c("GO:0044247","cellular polysaccharide catabolic process",0.099,1.8566,0.761,0.579,"oligosaccharide catabolism"),
                     c("GO:0006516","glycoprotein catabolic process",0.008,2.1149,0.796,0.327,"oligosaccharide catabolism"),
                     c("GO:0009311","oligosaccharide metabolic process",0.342,1.3576,0.863,0.508,"oligosaccharide catabolism"),
                     c("GO:0044265","cellular macromolecule catabolic process",1.268,1.4715,0.740,0.684,"oligosaccharide catabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
tiff( file="Results/RNAseq/revigo_treemap_REJ.tiff", res=300,width=2600, height=1900 ) # width and height are in inches
#pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Gene Ontology treemap genes up-regulated in REJ",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
