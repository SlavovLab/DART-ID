library(tidyverse)
source('Rscripts/lib.R')


## Fold change increase in IDs ------
#df <- fold.change.comp(list(RTLib=ev))

x <- logseq(1e-3, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'DART-ID', 'Percolator')
counter <- 0
for(i in x) {
  ratios <- c(
    1,
    sum(ev$pep_updated < i) /      sum(ev$PEP < i),
    sum(ev$pep_perc_updated < i) / sum(ev$PEP < i)
  )
  ident <- c(
    sum(ev$PEP < i) /              nrow(ev),
    sum(ev$pep_updated < i) /      nrow(ev),
    sum(ev$pep_perc_updated < i) / nrow(ev)
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=c('Spectra', 'DART-ID', 'Percolator'))

fold.change <- 
  ggplot(df) +
  geom_path(aes(x=x, y=(ratio-1)*100, color=Method), size=1) +
  geom_hline(aes(yintercept=1, color='Spectra'), size=1) +
  #geom_segment(aes(x=1e-2, xend=1e-2, y=-25, yend=180),
  #             color='black', linetype='dotted', size=0.5) +
  geom_vline(xintercept=1e-2, color='black', linetype='dotted', size=0.5) +
  scale_x_log10(limits=c(1e-3, 1), expand=c(0,0),
                breaks=c(1e-3, 1e-2, 1e-1, 1),
                labels=fancy_scientific) +
  scale_y_continuous(limits=c(-25, 170),
                     breaks=seq(0, 150, by=25),
                     expand=c(0,0)) +
  scale_color_manual(values=c(av[1], av[2], av[3])) +
  labs(x='Confidence Threshold', y='% Increase',
       title="Increase in confident PSMs",
       color=NULL) +
  theme_bert() + theme(
    plot.margin=unit(c(0.1,0.3,0.1,0.1), 'cm'),
    plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'), hjust=0.5),
    axis.title=element_text(size=10),
    legend.justification=c(0,1),
    legend.position=c(0.45, 1.06),
    legend.key.size=unit(0.35, 'cm'),
    legend.text=element_text(size=8),
    aspect.ratio=1
    #axis.line=element_line(color='black', size=0.5),
    #axis.text=element_text(size=12)
  )

# num.psms <- 
#   ggplot(df) +
#   geom_path(aes(x=x, y=ident, color=Method), size=1) +
#   #geom_vline(xintercept=1e-2, color='black', linetype='dashed', size=0.5) +
#   geom_segment(aes(x=1e-2, xend=1e-2, y=0, yend=0.68),
#                color='black', linetype='dotted', size=0.5) +
#   scale_x_log10(limits=c(1e-3, 1), 
#                 breaks=c(1e-3, 1e-2, 1e-1, 1), 
#                 labels=fancy_scientific,
#                 expand=c(0,0)) +
#   scale_y_continuous(limits=c(0,1),
#                      breaks=seq(0, 1, by=0.1),
#                      expand=c(0,0)) +
#   scale_color_manual(values=c(av[1], av[2], av[3])) +
#   labs(x='Confidence Threshold', y='Fraction of all PSMs', color=NULL,
#        title="Confident PSMs selected") +
#   theme_bert() + theme(
#     plot.margin=unit(c(0.1,0.3,0.1,0.1), 'cm'),
#     #axis.line=element_line(color='black', size=0.5),
#     legend.position=c(0.32, 0.87),
#     #legend.background = element_rect(fill='white'),
#     legend.key.size=unit(0.35, 'cm'),
#     legend.text=element_text(size=8),
#     plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'), hjust=0.5),
#     axis.title=element_text(size=10),
#     aspect.ratio=1
#   )

ggplot(df) +
  geom_path(aes(x=x, y=ident,     color='Spectra PEP'),    size=1) +
  geom_path(aes(x=x, y=ident_new, color='Spectra+RT PEP'), size=1) +
  #geom_vline(xintercept=1e-2, color='black', linetype='dashed', size=0.5) +
  geom_segment(aes(x=1e-2, xend=1e-2, y=0, yend=0.8),
               color='black', linetype='dotted', size=0.5) +
  scale_x_log10(limits=c(1e-3, 1), 
                breaks=c(1e-3, 1e-2, 1e-1, 1), 
                labels=fancy_scientific,
                expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  scale_color_manual(values=c(av[1], av[2])) +
  labs(x='PEP Threshold', y='Fraction', color=NULL,
       #title="Fraction of PSMs selected\n") +
       title="Data available for analysis\nafter filtering by confidence") +
  theme_bert() + theme(
    plot.margin=unit(c(0.1,0.5,0.1,0.1), 'cm'),
    #axis.line=element_line(color='black', size=0.5),
    legend.position=c(0.4, 0.91),
    #legend.background = element_rect(fill='white'),
    legend.key.size=unit(0.35, 'cm'),
    legend.text=element_text(size=8),
    plot.title=element_text(size=10, margin=margin(0,0,0,0,'cm'))
  )

#return(list(fold.change, num.psms))
