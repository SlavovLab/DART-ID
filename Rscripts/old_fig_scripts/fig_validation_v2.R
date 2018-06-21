library(tidyverse)
source('Rscripts/lib.R')
source('Rscripts/validate.lib.2.R')

## TMT Validation ------

# WARNING: this take a very, very long time!!!
cors <- validate.lib.2(ev, psm.freq.thresh=100)

# reorder levels
cors$Type = factor(cors$Type, levels(cors$Type)[c(3, 1, 2)])
cors$Type = plyr::revalue(cors$Type, c("New"="Spectra+PEP", "Original"="Spectra"))

fig2d <- 
  ggplot(cors) +
  #geom_line(stat='density', position='identity') +
  stat_density(aes(x=data, color=Type), adjust=3, 
               kernel='gaussian', geom='line', position='identity', size=1) +
  #geom_line(aes(x=))
  scale_color_manual(values=c(
    "Spectra+PEP"=av[2],
    "Spectra"=av[1],
    "Null"=av[5]
  )) +
  scale_x_continuous(expand=c(0, 0)) +
  labs(title=paste('Intra-Protein\nTMT Validation'),
       x='Correlations', y='Density', color=NULL, fill=NULL) +
  theme_bert() + theme(
    plot.margin=unit(c(0.1,0.5,0.1,0.1), 'cm'),
    axis.text=element_text(size=10),
    legend.key.height=unit(0.4, 'cm'),
    legend.position=c(0.3, 0.85)
  )

#return(fig2d)