library(tidyverse)

## Add percolator data --------

# ezconvert -v --config-file mq2pin -i ~/git/RTLib/Alignments/SQC_20180621_2/ev_updated.txt -o ~/git/RTLib/Alignments/SQC_20180621_2/pin.txt

# percolator -D 14 -U ~/git/RTLib/Alignments/SQC_20180621_2/pin.txt -m ~/git/RTLib/Alignments/SQC_20180621_2/pout.txt -Y

pout <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180724_3/pout.txt')
#pout <- read_tsv('/gd/bayesian_RT/Alignments/SQC2_20180815_1/pout.txt')

# line up IDs
pout <- pout %>% arrange(PSMId) %>%
  rename(perc_q_val=`q-value`)

#pout$PSMId %in% ev.f$id
ev <- cbind(ev, pout[match(ev$id, pout$PSMId),c('posterior_error_prob','perc_q_val')])

ev <- ev %>%
  mutate(pep_perc=posterior_error_prob,
         fdr_perc=perc_q_val) %>%
  mutate(pep_perc_updated=pep_perc)
ev$pep_perc_updated[is.na(ev$pep_perc)] = ev$PEP[is.na(ev$pep_perc)]