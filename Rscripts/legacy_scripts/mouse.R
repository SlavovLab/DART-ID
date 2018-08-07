library(tidyverse)
library(VIM)
source('Rscripts/lib.R')

# load data ---------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/mouse_scope_20180715_4/ev_updated.txt')

sort(unique(ev$`Raw file`))

dcols <- colnames(ev)[grepl('Reporter intensity corrected', colnames(ev))]
ev.f <- ev %>%
  # get UniProt accession ID
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))]) %>%
  # exclude CON, REV proteins
  filter(!grepl("CON__|REV__", Protein)) %>%
  # filter out sets with single cell equiv. lysates
  filter(!grepl('30[KL]|31A|32[HI]', `Raw file`)) %>%
  dplyr::select(c('Sequence', 'Modified sequence', 'Protein', 'Raw file', 'Retention time', 
                  'PEP', 'pep_new', 'pep_updated', 'qval', 'qval_updated', dcols))

ev.f <- ev.f %>%
  # remove empty and carrier channels
  dplyr::select(-grep('Reporter intensity corrected', colnames(ev.f))[c(8, 10)])


# further input processing ----------------------------------------------

# normalize data (by col, then by row)
dcols <- grep('Reporter intensity corrected', colnames(ev.f))
ev.f <- normalize_ri_data_table(ev.f, dcols, remove.empty.rows=F)

# remove rows without quantitation
# every exp, except for 43A-D, should have all cols filled out

empty_rows <- (apply(is.na(data.matrix(ev.f[,dcols])), 1, sum) > 0) &
  (!grepl('43[A-D]', ev.f$`Raw file`))
# for 43A-B, expect one empty col
empty_rows <- empty_rows | ((apply(is.na(data.matrix(ev.f[, dcols])), 1, sum) > 1) & (grepl('43[AB]', ev.f$`Raw file`)))
# for 43C-D, expect two empty cols
empty_rows <- empty_rows | ((apply(is.na(data.matrix(ev.f[, dcols])), 1, sum) > 2) & (grepl('43[CD]', ev.f$`Raw file`)))

sum(empty_rows)

ev.f <- ev.f[!empty_rows,]

# collapse quant by peptide and then by protein ---------------------------

dmat <- ev.f %>%
  filter(qval < 0.01) %>%
  group_by(`Raw file`, Protein, `Modified sequence`) %>%
  summarise_at(colnames(ev.f)[dcols], funs(mean)) %>%
  # collapse data by protein, by mean
  group_by(`Raw file`, Protein) %>%
  summarise_at(colnames(ev.f)[dcols], funs(mean))
  # throw out the protein names
  #dplyr::select(colnames(ev_a)[dcols])

dcols <- colnames(dmat)[grep('Reporter intensity corrected', colnames(dmat))]

# build channel-experiment pairs and spread them over the columns. 
# these represent single cells
dmat_a <- dmat %>% gather(channel, quant, dcols) %>%
  ungroup() %>%
  mutate(channel=paste(`Raw file`, channel, sep='_')) %>%
  dplyr::select(-c(`Raw file`)) %>%
  spread(channel, quant)
colnames(dmat_a)[-c(1)] <- seq(1:ncol(dmat_a)-1)

# remove proteins so we have at most 50% missing data
remove_prots <- apply(is.na(dmat_a), 1, sum) > (ncol(dmat_a) / 2)
dmat_a <- dmat_a[!remove_prots,]

# imputation --------------------------------------------------------------

dmat_a[,-c(1)] <- VIM::kNN(dmat_a[,-c(1)], imp_var=F)


# covariance --------------------------------------------------------------

heatmap(cor(dmat_a[,-c(1)]))

