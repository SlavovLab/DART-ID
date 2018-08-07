library(tidyverse)
library(reshape2)
library(pracma)
source('Rscripts/lib.R')

# increase in protein quant

# find out the indices of columns with reporter ion data
# should be 10 columns, but the code doesn't rely on this number
data.cols <- grep('Reporter intensity corrected', colnames(ev))
# ignore empty, carrier channels
# data.cols <- data.cols[c(5,6)]

# filter for sqc master sets only, and only keep correctly IDed proteins
ev.f <- ev %>%
  #filter(apply(ev[,data.cols]!=0, 1, sum) >= 8) %>%
  #filter(grepl('SQC', `Raw file`)) %>%
  filter(!grepl('SQC67[AB][16]|SQC67C1[3-9]|SQC67[CD]5|SQC68[DE]|IFN6[H-K]-Trg|SQC72D|SQC73[CD]|SQC74M|180416S_QC_SQC78A2',`Raw file`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!grepl('CON__', Proteins)) %>%
  # filter(apply(ev[,data.cols]!=0, 1, sum) == length(data.cols)) %>%
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  }))

ev.f <- ev.f %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
           seq(1, nrow(ev.f)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
           seq(1, nrow(ev.f)))[order(order(pep_updated))])

# only take rows w/ quantitation
# 8/10, because there are 2 empty channels
#ev.f <- ev.f[apply(ev.f[,data.cols] != 0, 1, sum) >= 8,]

# filter by FDR < 1%
fdr_thresh <- 0.01
ev.f.new  <- ev.f %>% filter(qval_updated < fdr_thresh)
ev.f.perc <- ev.f %>% filter(fdr_perc     < fdr_thresh)
ev.f      <- ev.f %>% filter(qval         < fdr_thresh)

experiments <- sort(unique(ev.f$`Raw file`))
#prots <- unique(ev.f.new$Protein)
#prots <- unique(ev$Protein)
prots <- unique(ev$`Modified sequence`)

dmat <- zeros(length(prots), length(experiments))
dmat_new <- zeros(length(prots), length(experiments))
dmat_perc <- zeros(length(prots), length(experiments))

# Proteins
# for(i in 1:length(experiments)) {
#   cat('\r', i, '/', length(experiments), '         ')
#   flush.console()
#
#   ev.a <- ev.f %>% filter(`Raw file`==experiments[i])
#   dmat     [prots %in% ev.a$Protein,i] <- 1
#   ev.a <- ev.f.new %>% filter(`Raw file`==experiments[i])
#   dmat_new [prots %in% ev.a$Protein,i] <- 1
#   ev.a <- ev.f.perc %>% filter(`Raw file`==experiments[i])
#   dmat_perc[prots %in% ev.a$Protein,i] <- 1
# }

# Peptides
for(i in 1:length(experiments)) {
  cat('\r', i, '/', length(experiments), '         ')
  flush.console()
  
  ev.a <- ev.f %>% filter(`Raw file`==experiments[i])
  dmat     [prots %in% ev.a$`Modified sequence`,i] <- 1
  ev.a <- ev.f.new %>% filter(`Raw file`==experiments[i])
  dmat_new [prots %in% ev.a$`Modified sequence`,i] <- 1
  ev.a <- ev.f.perc %>% filter(`Raw file`==experiments[i])
  dmat_perc[prots %in% ev.a$`Modified sequence`,i] <- 1
}

# dmat_c
# - 0 if not quantified
# - 1 if quantified with FDR < 1%
# - 2 if quantified with FDR > 1% and FDR_new < 1%
# - 3 if no longer quantified with FDR < 1% and FDR_new > 1%
dmat_c <- dmat
dmat_c[!dmat & dmat_new] <- 2
dmat_c[dmat & !dmat_new] <- 3

# reorder the peptides/proteins so that the 1s are stacked on the bottom
#odr <- rev(order(apply(dmat, 1, sum) + apply(dmat_new, 1, sum)))
odr <- rev(order(apply(dmat, 1, sum)))
dmat_c <- dmat_c[odr,]

# remove rows with all 0s
dmat_c <- dmat_c[apply(dmat_c, 1, sum) != 0,]


