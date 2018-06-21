validate.lib <- function(
  ev.in, 
  ev.type='STAN', 
  cell.types=c('J', 'H'),     # subset of cell types to compare
  exps=NULL,                  # subset of experiments to use
  remove.REV=TRUE,            # remove REV (reverse) matches
  remove.CON=TRUE,            # remove CON (contaminants)
  unique.only=TRUE,           # only take unique peptides
  pep.freq.thresh=10,         # minimum number of peptides per protein
  cor.size.thresh=10,         # minimum number of observations/PSMs per protein
                              # to build the correlation matrix
  melt.output=FALSE           # melt the data frame?
) {

  library(tidyverse)
  library(reshape2)
  library(stringr)
  source('lib.R')
  
  # load evidence file
  ev <- parse.ev.adj(ev.in, ev.type)
  
  # assign protein IDs
  ev$Protein_ID <- as.numeric(as.factor(ev$Proteins))
  
  # clean up raw file names and extract an experiment ID from it
  # need the experiment ID (19A, 30B, etc.), to match with sample metadata from the
  # experiment description excel sheet/.csv
  ev$file <- clean.file.name(ev$`Raw file`)
  ev$exp <- str_extract(ev$file, '[1-9][0-9][A-Z](\\d)?')
  
  # load excel description
  desc <- process.desc()
  desc.types <- dcast(desc, Exp~Type, value.var='Type')
  desc.types$Exp <- as.character(desc.types$Exp)
  
  # take subset of experiments
  # choose experiments with at least one of the specified cell types present
  # cell.types defaults to c('J', 'H'),
  # for Jurkat and Hek cells
  exps.cell.type <- desc.types[apply(desc.types[,cell.types], 1, sum) > 1,'Exp']
  # if specified, intersect with list of experiments given by user
  if(is.null(exps)) {
    exps <- exps.cell.type
  } else {
    exps <- intersect(exps, exps.cell.type)
  }
  # make sure we actually have experiments to work with
  if(length(exps) < 1) {
    stop(paste0('No experiments given, or meet cell type requirements\n',
                'Experiments must have cell types: ', paste(cell.types, collapse=', ')))
  }
  ev <- ev[ev$exp %in% exps,]
  
  # find out the indices of columns with reporter ion data
  # should be 10 columns, but the code doesn't rely on this number
  data.cols <- grep('intensity', colnames(ev))
  
  # remove peptides without a protein
  # this can happen, especially when searching with a restricted FASTA file
  ev <- ev[!is.na(ev$Proteins),]
  # remove non-unique peptides
  # aka, peptides with only one entry in the proteins column
  if(unique.only) {
    ev <- ev[!grepl(';',ev$Proteins),]
  } else {
    # otherwise, default to razor protein
    ev$Proteins <- ev$`Leading razor protein`
  }
  
  # remove REV proteins
  # these are MaxQuant decoys that scored higher than the actual
  # peptide of interest
  if(remove.REV) {
    ev <- ev[!grepl('REV*',ev$Proteins),]
  }
  # remove CON proteins?
  # common contaminants, like Keratin, Albumin, etc.
  if(remove.CON) {
    ev <- ev[!grepl('CON*',ev$Proteins),]
  }
  
  # set zeroes to NA
  ev[,data.cols][ev[,data.cols]==0] = NA
  # remove rows without quantitation
  ev <- ev[!(apply(ev[,data.cols], MARGIN=1, FUN=sum, na.rm=TRUE)==0),]
  
  ## remove carrier and empty channels
  # these channel indices vary between experiments, 
  # so we're gonna have to loop thru these and check for each one
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
  
  ## only use proteins with many peptides
  
  # count # of peptides per protein in our entire dataset
  prot.map <- as.data.frame(table(ev$Proteins))
  prot.map <- prot.map[order(prot.map$Freq, decreasing=TRUE),]
  
  # only use proteins with more than [pep.freq.thresh] peptides attached to it
  # pep.freq.thresh defaults to 10
  prots <- as.character(prot.map$Var1[prot.map$Freq >= pep.freq.thresh])
  
  ## Correlate peptides within each protein
  #
  # Iterate over proteins and do pairwise correlation between 
  # quantitation vectors of peptides.
  # Do this for certain subsets of the protein data, 
  # in order to compare them after the fact
  
  print('Building correlation vectors...\n')
  cors <- as.data.frame(t(sapply(prots, FUN=function(prot) {
    cat('\r', match(prot, prots), '/', length(prots), 
        '                           ')
    flush.console()
    
    # get PSMs for this protein
    
    # ORIGINAL 
    # -- originally good IDs (PEP < 0.05), ignoring result of RTLib method
    prot.data <- ev[ev$Proteins==prot & ev$PEP < 0.05,data.cols]
    
    # NEW
    # -- originally bad IDs, upgraded to good ID (PEP > 0.05 & PEP.new < 0.05)
    # -- set of NEW and ORIGINAL should be disjoint
    prot.data.new <- ev[ev$Proteins==prot & ev$PEP > 0.05 & ev$PEP.new < 0.05,data.cols]
    
    # TOTAL
    # -- good IDs after bayesian update, regardless of if it was upgraded beyond threshold or not
    # -- **NOT** the union of ORIGINAL and NEW (but does contain all of NEW), 
    # -- as some PSMs will be downgraded below PEP threshold via. the RTLib method
    # --
    # -- reminder: PEP.updated is PEP superceded by PEP.new, so that there are no NAs
    # -- that is, a PSM that was updated by the bayesian update will have its PEP replaced by PEP.new
    prot.data.tot <- ev[ev$Proteins==prot & ev$PEP.updated < 0.05,data.cols]
    
    # init output
    prot.cor <- NULL
    prot.cor.new <- NULL
    prot.cor.tot <- NULL
    
    # lets say... must have more than [cor.size.thresh] observations to have a meaninful correlation matrix
    # not enough data can push extremely low/high correlations, especially when so many
    # of our columns are NA anyways
    # cor.size.thresh defaults to 10
    #
    # the pairwise.complete.obs method of the cor() function can exacerbate 
    # this effect if there aren't enough observations
    
    if(nrow(prot.data) > cor.size.thresh) {
      # pairwise complete to account to missing/NA values
      # results in matrix that is not positive semi-definite, and rows that may contain NA
      prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
      # set diagonal to NA
      diag(prot.cor) <- NA
    }
    # do the same above, but for prot.data.new, prot.data.tot
    if(nrow(prot.data.new) > cor.size.thresh) {
      prot.cor.new <- cor(t(prot.data.new), use='pairwise.complete.obs', method='pearson')
      diag(prot.cor.new) <- NA
    }
    if(nrow(prot.data.tot) > cor.size.thresh) {
      prot.cor.tot <- cor(t(prot.data.tot), use='pairwise.complete.obs', method='pearson')
      diag(prot.cor.tot) <- NA
    }
    
    # take the median of the correlation matrix and pass it up
    # sometimes either prot.data, new, or tot will have no rows
    # in that case just output NA and have the density() function
    # or whatever is evaluating this later exclude it from calculations
    c(ifelse(nrow(prot.data) > cor.size.thresh, median(prot.cor, na.rm=TRUE), NA),
      ifelse(nrow(prot.data.new) > cor.size.thresh, median(prot.cor.new, na.rm=TRUE), NA),
      ifelse(nrow(prot.data.tot) > cor.size.thresh, median(prot.cor.tot, na.rm=TRUE), NA))
  })))
  
  # clean up the resultant data frame for cors
  cors$Protein <- rownames(cors)
  rownames(cors) <- NULL
  colnames(cors) <- c('Original', 'New', 'Total', 'Protein')
  cors <- cors[,c(4, 1, 2, 3)]
  
  ## Compute null distribution for peptide/protein correlations
  
  cors$Null <- NA
  # set seed for consistency
  #set.seed(1)
  # take n at a time
  # n = median # of observations for this set of proteins
  # this should prevent too many/too little observations from affecting
  # the null distribution vs. the correct ones
  n = median(prot.map$Freq[prot.map$Var1 %in% prots])
  for(i in 1:nrow(cors)) {
    ## pairwise correlation between peptides from random proteins
    
    # only use proteins that were used in the previous correlations
    #prot.inds <- which(ev$Proteins %in% prots)
    
    # randomly sample peptides from this list, and then compute the
    # pairwise correlation between them
    #
    # ideally this will be uniform - but will probably never be as
    # that would assume that all proteins are regulated independently
    # there are definitely proteins whose levels will correlate closely with others
    prot.inds <- which(ev$PEP.updated > 0.5)
    prot.data <- ev[sample(prot.inds, size=n), data.cols]
    
    # nevermind, just randomly match PSMs, 
    # regardless of parent protein or whether or not it was used in the previous set
    #prot.data <- ev[sample(1:nrow(ev), size=n), data.cols]
    
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    
    # set diagonal to NA
    diag(prot.cor) <- NA
    # take the median of the correlation matrix and pass it up
    cors$Null[i] <- median(prot.cor, na.rm=TRUE)
  }
  
  if(melt.output) {
    # melt into a form that ggplot2 likes
    cors.m <- melt(cors, id.var='Protein', variable.name='type', value.name='Correlation')
    return(cors.m)
  }
  
  return(cors)
  
  
  #cors.m %>%
  #  ggplot(aes(x=Correlation, color=type)) +
  #  geom_density()
}