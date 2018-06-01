library(tidyverse)
source('lib.R')
source('align.rt.R')

ev <- read_tsv("/gd/Slavov_Lab/SingleCell_Data/FP17/evidence.txt")
ev <- ev %>%
  rename(`Sequence ID`=`Peptide ID`) %>%
  rename(`Peptide ID`=`Mod. peptide ID`) # alias the modified peptide ID as the peptide ID

#align.rt(ev, pars.out='params.FP17.RData', figs.out='./FP17_align')
align.rt(ev, pars.out='params.FP17.RData', figs.out='./FP17_align', STAN.iters=1e5, rt.distortion=1, prior.iters=10, pep_thresh=0.5, rtl_filter=1)

load('params.FP17.RData')

adjust.pep(ev, pars, pep_thresh=0.5, rtl_filter=1, n_exp=3, out.path='FP17.ev.adj.txt')

ev.a <- data.frame()
for(exp in 1:num_experiments) {
  exp_indices=ev.f$exp_id==exp
  ev.b <- subset(ev.f, exp_indices)
  
  
  predicted <- (muij_fit[muij_map])[exp_indices[muij_map]]
  observed <- retention_times[exp_indices[muij_map]]
  residual <- observed - predicted
  
  inds <- as.numeric(which(log(abs(residual)) > 1))
  ev.c <- ev[match(ev.b[inds,]$`Peptide ID`, ev$`Peptide ID`),]
  
  ev.a <- rbind(ev.a, data.frame(Protein=as.character(ev.c$`Leading razor protein`), residual=as.numeric(log(abs(residual[inds]))), pep=as.numeric(ev.c$PEP)))
}


# python
#ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180506_mu_norm_1/ev_adjusted.txt")
ev <- read_tsv("~/git/RTLib/Alignments/FP17_20180507_2/ev_adjusted.txt")
ev <- ev %>% rename(PEP.new=pep_new)
# R
ev <- read_tsv("FP17.ev.adj.txt")
source('validate.lib.R')
source('validate.lib.2.R')

cors <- validate.lib(ev)
cors <- validate.lib.2(ev)

ggplot(cors, aes(x=data, color=Type)) +
  #geom_histogram(aes(y=..density.., color=NULL, fill=Type), position='identity', alpha=0.5, bins=30) +
  geom_line(stat='density', position='identity') +
  scale_color_manual(values=c(
    "New"="blue",
    "Original"="black",
    "Null"="red"
  )) +
  scale_fill_manual(values=c(
    "New"="blue",
    "Original"="black",
    "Null"="red"
  ), guide=FALSE) +
  labs(title=paste('FP17 Validation'),
       x='Correlations', y='Density', color=NULL, fill=NULL) +
  theme_bert()
