library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

## ---------

# melt into a frame that can be plotted by ggplot
dmat_cc <- melt(dmat_c)
dmat_cc$value <- as.factor(dmat_cc$value)

## ---------

ggplot(dmat_cc, aes(x=Var2, y=Var1, fill=value)) +
geom_raster(interpolate=F) +
scale_x_continuous(limits=c(0, length(experiments)), 
                   breaks=seq(0, length(experiments), by=20), 
                   expand=c(0,0)) +
scale_y_continuous(limits=c(0, nrow(dmat_c)), 
                   breaks=seq(0, length(prots), by=500),
                   expand=c(0,0)) +
scale_fill_manual(values=c("#FFFFFF", "#000000", "#FF0000", "#0000FF"), 
                  labels=c("Not Quantified", 
                           "Quantified (Spectra PEP < 0.01)", 
                           "Quantified (DART-ID PEP < 0.01)",
                           'Quantified (not anymore)')) +
labs(x="Single Cell Experiments", y="Peptide Sequence", fill=NULL,
     title="Peptide Coverage Increase") +
guides(fill=guide_legend(override.aes=list(color="black"), nrow=3)) +
theme_bert() + theme(
  plot.margin=unit(c(0.1,0.3,0.1,0.1), 'cm'),
  plot.title=element_text(size=10, margin=margin(0,0,0.05,0,'cm'), hjust=0.5),
  axis.title=element_text(size=10),
  #axis.line=element_line(color='black', size=1),
  legend.position='bottom',
  legend.key=element_rect(color='black'),
  legend.text=element_text(size=8),
  legend.key.height=unit(0.3, 'cm'),
  legend.key.width=unit(0.3, 'cm'),
  legend.margin=margin(0,0,0,0, unit='cm'),
  legend.box.spacing=unit(c(0.1,0,0,0), 'cm')
)
