library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')

## load data -------

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
  # only select SQC master sets
  filter(grepl("SQC", `Raw file`)) %>%
  # filter out non J/U sets (QC40E,F, QC42A)
  filter(!grepl('SQC9', `Raw file`)) %>%
  dplyr::select(c('Sequence', 'Modified sequence', 'Protein', 'Raw file', 'Retention time', 
           'PEP', 'pep_new', 'pep_updated', 'qval', 'qval_updated', dcols)) %>% 
  # remove empty channels
  dplyr::select(-grep('Reporter intensity corrected', colnames(ev.f))[c(3, 4)])
#  filter(qval < 0.01) %>%

# only take rows w/ quantitation
ev.f <- ev.f[apply(ev.f[,grep('Reporter', colnames(ev.f))] == 0, 1, sum) == 0,]

# normalize data (by col, then by row)
dcols <- grep('Reporter intensity corrected', colnames(ev.f))
ev.f <- normalize_ri_data_table(ev.f, dcols)

## -------

conf_thresh <- 0.001

# prot_quants <- ev.f %>% 
#   #filter(!is.na(pep_new)) %>%
#   group_by(Protein, `Raw file`) %>%
#   summarise(confident=any(qval < 0.01),
#             new_confident=any(qval > 0.01 & qval_updated < 0.01)) %>%
#   group_by(Protein) %>%
#   summarise(num_confident=sum(confident),
#             num_new_confident=sum(new_confident)) %>%
#   mutate_at(c('num_confident', 'num_new_confident'), funs(ifelse(is.na(.), 0, .))) %>%
#   filter(num_confident != 0 | num_new_confident != 0)

# old_prots <- prot_quants %>% filter(num_confident > 0) %>% pull(Protein)
# new_prots <- prot_quants %>% filter(num_confident == 0) %>% pull(Protein)

old_prots <- sort(unique(ev.f %>% filter(qval < conf_thresh) %>% pull(Protein)))
new_prots <- sort(unique(ev.f %>% filter(qval_updated < conf_thresh) %>% pull(Protein)))
#old_prots <- sort(unique(ev.f %>% filter(PEP < 0.01) %>% pull(Protein)))
#new_prots <- sort(unique(ev.f %>% filter(pep_updated < 0.01) %>% pull(Protein)))
new_prots <- new_prots[!new_prots %in% old_prots]

## filter, normalize, collapse, cluster data -------

ev_a <- ev.f %>%
  # only select from protein list
  filter(Protein %in% old_prots) %>%
  # filter at 1% FDR
  filter(qval < conf_thresh)

ev_b <- ev.f %>%
  filter(Protein %in% new_prots) %>%
  filter(qval_updated < conf_thresh)

dmat_a <- ev_a %>%
  # collapse data by sequence, by mean
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_a)[dcols], funs(mean)) %>%
  # collapse data by protein, by mean
  group_by(Protein) %>%
  summarise_at(colnames(ev_a)[dcols], funs(mean)) %>%
  # throw out the protein names
  dplyr::select(colnames(ev_a)[dcols])

dmat_b <- ev_b %>%
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_b)[dcols], funs(mean)) %>%
  group_by(Protein) %>%
  summarise_at(colnames(ev_b)[dcols], funs(mean)) %>%
  dplyr::select(colnames(ev_b)[dcols])

# cast to matrix, and manually cluster (J, then U)
dmat_a <- data.matrix(dmat_a)
dmat_b <- data.matrix(dmat_b)
dmat_a <- dmat_a[,c(1, 3, 5, 7, 2, 4, 6, 8)]
dmat_b <- dmat_b[,c(1, 3, 5, 7, 2, 4, 6, 8)]

## get correlation matrices ------

cor_type_a <- cor(dmat_a)
cor_type_b <- cor(dmat_b)

j_channels <- c(2, 3, 4)
u_channels <- c(6, 7, 8)

ratio_mat_a <- zeros(nrow(dmat_a), length(j_channels) * length(u_channels))
ratio_mat_b <- zeros(nrow(dmat_b), length(j_channels) * length(u_channels))
for(j in 1:length(j_channels)) {
  for(u in 1:length(u_channels)) {
    ratio_mat_a[,((j-1)*length(j_channels))+u] <- 
      dmat_a[,j_channels[j]] / dmat_a[,u_channels[u]]
    ratio_mat_b[,((j-1)*length(j_channels))+u] <- 
      dmat_b[,j_channels[j]] / dmat_b[,u_channels[u]]
  }
}

cor_ratio_a <- cor(log2(ratio_mat_a))
cor_ratio_b <- cor(log2(ratio_mat_b))
# cor_ratio_a <- cor(ratio_mat_a)
# cor_ratio_b <- cor(ratio_mat_b)

# ## cell-type correlations --------
# 
# colfunc <- colorRampPalette(c('red', 'white', 'blue'))
# image(cor(dmat),col=colfunc(20))
# 
# #image.plot(cor(dmat), breaks=seq(-1, 1, length.out=21), col=colfunc(20))
# 
# ## single-cell ratios --------
# j_channels <- c(2, 3, 4)
# u_channels <- c(6, 7, 8)
# 
# ratio_mat <- zeros(nrow(dmat), length(j_channels) * length(u_channels))
# for(j in 1:length(j_channels)) {
#   for(u in 1:length(u_channels)) {
#     ratio_mat[,((j-1)*length(j_channels))+u] <- dmat[,j_channels[j]] / dmat[,u_channels[u]]
#   }
# }
# 
# image(cor(ratio_mat), col=colfunc(100))

