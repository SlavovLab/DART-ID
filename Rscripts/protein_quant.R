library(tidyverse)
library(reshape2)
source('Rscripts/lib.R')

# increase in protein quant

# find out the indices of columns with reporter ion data
# should be 10 columns, but the code doesn't rely on this number
data.cols <- grep('Reporter intensity corrected', colnames(ev))
# ignore empty, carrier channels
# data.cols <- data.cols[c(5,6)]

# extract uniprot ID
ev$Protein <- sapply(strsplit(ev$`Leading razor protein`, "\\|"), function(p) {
  if(length(unlist(p)) == 1) return(p[1])
  else if(length(unlist(p)) == 3) return(p[2])
  else return(p[1])
})

# filter for sqc master sets only, and only keep correctly IDed proteins
ev.f <- ev %>%
  # filter(apply(ev[,data.cols]!=0, 1, sum) == length(data.cols)) %>%
  filter(apply(ev[,data.cols]!=0, 1, sum) >= 8) %>%
  filter(!grepl("IFN|FP18", `Raw file`)) %>%
  filter(!grepl("REV__", `Leading razor protein`)) %>%
  filter(!grepl("CON__", Proteins))

# filter by PEP < 0.01 (TODO - do this by FDR)
pep_thresh <- 0.01
ev.f.new  <- ev.f %>% filter(pep_updated < pep_thresh)
ev.f.perc <- ev.f %>% filter(pep_perc    < pep_thresh)
ev.f      <- ev.f %>% filter(PEP         < pep_thresh)

experiments <- sort(unique(ev.f$`Raw file`))
#prots <- unique(ev.f.new$Protein)
#prots <- unique(ev$Protein)
prots <- unique(ev$Sequence)

dmat <- zeros(length(prots), length(experiments))
dmat_new <- zeros(length(prots), length(experiments))
dmat_perc <- zeros(length(prots), length(experiments))

# Proteins
# for(i in 1:length(experiments)) {
#   ev.a <- ev.f %>% filter(`Raw file`==experiments[i])
#   dmat     [prots %in% ev.a$Protein,i] <- 1
#   ev.a <- ev.f.new %>% filter(`Raw file`==experiments[i])
#   dmat_new [prots %in% ev.a$Protein,i] <- 1
#   ev.a <- ev.f.perc %>% filter(`Raw file`==experiments[i])
#   dmat_perc[prots %in% ev.a$Protein,i] <- 1
# }

# Peptides
for(i in 1:length(experiments)) {
  ev.a <- ev.f %>% filter(`Raw file`==experiments[i])
  dmat     [prots %in% ev.a$Sequence,i] <- 1
  ev.a <- ev.f.new %>% filter(`Raw file`==experiments[i])
  dmat_new [prots %in% ev.a$Sequence,i] <- 1
  ev.a <- ev.f.perc %>% filter(`Raw file`==experiments[i])
  dmat_perc[prots %in% ev.a$Sequence,i] <- 1
}

# dmat_c
# - 0 if not quantified
# - 1 if quantified with PEP < 0.01
# - 2 if quantified with PEP > 0.01 and PEP_new < 0.01
# - 3 if no longer quantified with PEP < 0.01 and PEP_new > 0.01
dmat_c <- dmat
dmat_c[!dmat & dmat_new] <- 2

# reorder the peptides/proteins so that the 1s are stacked on the bottom
#odr <- rev(order(apply(dmat, 1, sum) + apply(dmat_new, 1, sum)))
odr <- rev(order(apply(dmat, 1, sum)))
dmat_c <- dmat_c[odr,]

# remove rows with all 0s
dmat_c <- dmat_c[apply(dmat_c, 1, sum) != 0,]

# melt into a frame that can be plotted by ggplot
dmat_cc <- melt(dmat_c)
dmat_cc$value <- as.factor(dmat_cc$value)