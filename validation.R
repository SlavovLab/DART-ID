library(pracma)
library(tidyverse)
library(reshape2)
library(stringr)
require(RODBC)
source('lib.R')
ev <- parse.ev.adj('dat/ev.adj.Fit2.txt')
# assign protein IDs
ev$Protein_ID <- as.numeric(as.factor(ev$Proteins))

ev$file <- ev$`Raw file`
ev$file <- regexprep(ev$file, '#', '')
ev$file <- regexprep(ev$file, '.*_NC_', '')
ev$file <- regexprep(ev$file, '.raw', '')
ev$file <- regexprep(ev$file, 'set', '')
ev$file <- regexprep(ev$file, '[\\-+=%]', '_')
ev$exp <- str_extract(ev$file, '[1-9][0-9][A-Z](\\d)?')

# load excel description
desc <- process.desc()
desc.types <- dcast(desc, Exp~Type, value.var='Type')
desc.types$Exp <- as.character(desc.types$Exp)

# choose experiments with at least one J and one H channel
exps <- desc.types[desc.types$J > 1 & desc.types$H > 1,'Exp']
ev <- ev[ev$exp %in% exps,]

data.cols <- grep('intensity', colnames(ev))
# ignore empty and carrier channels - we won't use them
#data.cols <- data.cols[-c(8, 10)]

# remove peptides w/o a protein
ev <- ev[!is.na(ev$Proteins),]
# remove non-unique peptides
ev <- ev[!grepl(';',ev$Proteins),]

# remove REV proteins
ev <- ev[!grepl('REV*',ev$Proteins),]
# remove CON proteins?
ev <- ev[!grepl('CON*',ev$Proteins),]
# remove PEP > 0.05
ev <- ev[ev$PEP > 0.05,]

# take subset of experiments
ev <- ev[grepl('30[A-J]|29[A-C]|26[A-E]|25[A-C]|24[A-C]', ev$exp),]

# set 0 to NA
ev[,data.cols][ev[,data.cols]==0] = NA

# remove rows without quantitation
#ev <- data[!(apply(data, 1, FUN=sum, na.rm=TRUE)==0),]
ev <- ev[!(apply(ev[,data.cols], MARGIN=1, FUN=sum, na.rm=TRUE)==0),]

## remove carrier and empty channels
## these channels vary between experiments, so we're gonna have to
## loop thru these
for(i in exps) {
  if(sum(ev$exp==i) <= 0) next
  inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | is.na(desc$Type)),'ch']
  # make relative to the ev data frame
  inds.remove <- inds.remove + data.cols[1] - 1
  # set to NA
  ev[ev$exp==i,inds.remove] <- NA
}

## normalize data

# first normalize by column, by median
# this assumes the same amount of protein per channel
# (not by cell, but in aggregate - this erases any ionization or 
# collision cell or any other quantitative biases)
ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], MARGIN=2, FUN=median, na.rm=TRUE)

# now normalize across rows, to get the difference between the channels
# AKA the difference between different cells/conditions
ev[,data.cols] <- ev[,data.cols] / t(apply(ev[,data.cols], MARGIN=1, FUN=median, na.rm=TRUE))

# take median of J cells and median of H cells
for(i in exps) {
  cat('\r', match(i, exps), '/', length(exps), '  ')
  flush.console()
  #ev[ev$exp==i]
  inds <- ev$exp==i
  if(sum(inds) <= 0) next
  # only take single cell channels
  j.inds <- desc[desc$Exp==i & desc$Type=='J' & 
                   !is.na(desc$Type) & desc$Quantity < 10,'ch'] + data.cols[1] - 1
  h.inds <- desc[desc$Exp==i & desc$Type=='H' & 
                   !is.na(desc$Type) & desc$Quantity < 10,'ch'] + data.cols[1] - 1
  ev[inds,'J'] <- apply(ev[inds,j.inds], MARGIN=1, FUN=median, na.rm=TRUE)
  ev[inds,'H'] <- apply(ev[inds,h.inds], MARGIN=1, FUN=median, na.rm=TRUE)
  ev[inds, 'JHvar'] <- apply(ev[inds,data.cols], MARGIN=1, FUN=var, na.rm=TRUE)
}
# remove rows that don't have J and H to compare against
ev <- ev[apply(is.na(ev[,c('J','H')]), MARGIN=1, FUN=sum)==0,]

# count # of peptides per protein in our entire dataset
prot.map <- as.data.frame(table(ev$Proteins))
prot.map <- prot.map[order(prot.map$Freq, decreasing=TRUE),]

# only use proteins with more than 10 peptides attached to it
prots <- as.character(prot.map$Var1[prot.map$Freq >= 20])

cors <- sapply(prots, FUN=function(prot) {
  cat('\r', match(prot, prots), '/', length(prots), '                           ')
  flush.console()
  
  # get peptide data for this protein
  #prot.data <- ev[ev$Proteins==prot,c('J','H')]
  prot.data <- ev[ev$Proteins==prot,data.cols]
  # pairwise complete to account to missing/NA values
  # results in matrix that is not positive semi-definite, and rows that may contain NA
  prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
  #prot.cor <- as.numeric(cor(prot.data[,1], prot.data[,2], use='complete.obs', method='pearson'))
  
  # set diagonal to NA
  #diag(prot.cor) <- NA
  
  # remove rows/cols with all NAs
  #prot.inds.remove <- which(apply(prot.cor, MARGIN=1, FUN=sum, na.rm=TRUE)==0)
  # remove rows that contain NA after this (unlucky pairwise comparisons), 
  # not counting the rows that were all NAs in the first place
  # keep track for indexing purposes
  #prot.inds.remove <- sort(c(prot.inds.remove, 
  #  which(apply(is.na(prot.cor), MARGIN=1, FUN=sum) > length(prot.inds.remove))))
  
  # reform the correlation matrix
  #prot.cor <- prot.cor[-prot.inds.remove, -prot.inds.remove]
  # only take the upper diagonal, and form into vector
  #prot.cor <- prot.cor[upper.tri(prot.cor)]
  median(prot.cor, na.rm=TRUE)
  
  #prot.cor
})

# get null distribution for correlation
# pairwise correlation between peptides from random proteins
null.cors <- vector(length=length(cors))
par(mfrow=c(2,2))
for(i in 1:4) {
  # take 50 at a time
  # randomize this quantity????
  n = 50
  for(i in 1:length(cors)) {
    prot.data <- ev[sample.int(nrow(ev), size=n), data.cols]
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    # set diagonal to NA
    diag(prot.cor) <- NA
    null.cors[i] <- median(prot.cor, na.rm=TRUE)
  }
  
  #plot(density(null.cors), col='red')
  #lines(density(cors), col='blue')
  plot(density(cors), col='blue')
  lines(density(null.cors), col='red')
  #lines(c(0,0),c(0,1e5),type='l')
}
par(mfrow=c(1,1))





ev.a <- ev[ev$Proteins=='sp|P13639|EF2_HUMAN',]
ev.cor <- cor(t(ev.a[,data.cols]), use='pairwise.complete.obs')
#ev.cor <- cor(t(ev.a[,data.cols]), use='everything')
# remove rows/cols with all NAs
inds.remove <- which(apply(ev.cor, MARGIN=1, FUN=sum, na.rm=TRUE)==0)

# remove rows that contain NA after this (unlucky pairwise comparisons), 
# not counting the rows that were all NAs in the first place
# keep track for indexing purposes
inds.remove <- sort(c(inds.remove, 
                     which(apply(is.na(ev.cor), MARGIN=1, FUN=sum) > length(inds.remove))))

ev.cor.dist <- as.dist((1-ev.cor[-inds.remove,-inds.remove])^2)
ev.hclust <- hclust(ev.cor.dist)

image(as.matrix(ev.cor.dist)[ev.hclust$order,ev.hclust$order])
  