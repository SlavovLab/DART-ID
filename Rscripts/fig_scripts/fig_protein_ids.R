library(tidyverse)
source('Rscripts/lib.R')


df <- data.frame(x=factor(c("Spectra", "DART-ID", "Percolator"), 
                          levels=c("Spectra", "DART-ID", "Percolator")),
                 y=c(sum(apply(dmat, 1, sum)) / length(experiments), 
                     sum(apply(dmat_new, 1, sum)) / length(experiments),
                     sum(apply(dmat_perc, 1, sum)) / length(experiments)))
fig2f <- 
  ggplot(df,aes(x=x, y=y, fill=x)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c(av[1], av[2], av[3]), guide=F) +
  scale_y_continuous(limits=c(0, 1200), breaks=seq(0, 1200, by=200), expand=c(0,0)) +
  labs(x=NULL, y="Peptides Quantified per Experiment", fill=NULL, title="\n") +
  theme_bert() + theme(
    plot.margin=unit(c(0.1,0.1,0.1,0.1), 'cm'),
    #axis.line = element_line(color='black', size=0.75),
    axis.title=element_text(size=10),
    axis.text.x = element_text(angle=45, hjust=1),
    plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'))
  )

# violin plot?
df <- data.frame(
  x=factor(rep(c("Spectra", "DART-ID", "Percolator"), each=ncol(dmat)),
           levels=c("Spectra", "DART-ID", "Percolator")),
  y=as.numeric(c(apply(dmat, 2, sum), apply(dmat_new, 2, sum), apply(dmat_perc, 2, sum)))
)
#fig2f <- 
# ggplot(df, aes(x, y, group=x, fill=x)) +
#   #geom_violin() +
#   geom_boxplot(outlier.color='white') +
#   scale_y_continuous(limits=c(0, 2000), breaks=seq(0, 2000, by=500), expand=c(0,0)) +
#   scale_fill_manual(values=c(av[1], av[2], av[3]), guide=F) +
#   labs(x=NULL, y="Peptides Quantified per Experiment", fill=NULL, title="\n") +
#   theme_bert() + theme(
#     plot.margin=unit(c(0.1,0.1,0.1,0.1), 'cm'),
#     #axis.line = element_line(color='black', size=0.75),
#     axis.title=element_text(size=10),
#     axis.text.x = element_text(angle=45, hjust=1),
#     plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'))
#   )