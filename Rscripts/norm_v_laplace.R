library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
ev_n <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_2/ev_updated.txt")

ev.f <- ev %>%
  filter(!is.na(pep_new))

ev_n.f <- ev_n %>%
  filter(!is.na(pep_new))

res <- abs(ev.f$`Retention time` - ev.f$muij)
res_n <- abs(ev_n.f$`Retention time` - ev_n.f$muij)

plot(density(res, n=1e4), xlim=c(0, 2))
lines(density(res_n, n=1e4), col='blue')
