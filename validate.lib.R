validate.lib <- function(
  ev, 
  ev.type='STAN',
  #exps=NULL,                 # subset of experiments to use
  raw.files=NULL,             # subset of raw files to use
  remove.REV=TRUE,            # remove REV (reverse) matches
  remove.CON=TRUE,            # remove CON (contaminants)
  remove.keratins=TRUE,       # remove keratins (using Keratins.txt file)
  sparse.filter=0.5,          # filter out sparse observations. 
                              # this is the fraction of columns in 
                              # the TMT data that is allowed to be NA
  unique.only=FALSE,          # only take unique peptides
  norm.cols=TRUE,             # normalize by columns 
                              # assume same amt of protein per channel
  norm.rows=TRUE,             # normalize by rows 
                              # assume same amt of protein per observation
  psm.freq.thresh=10,         # minimum number of PSMs per protein
  alpha=0.05,                 # significance threshold that separates
                              # Original and New sets (PEP < alpha)
  cor.size.thresh=5           # minimum number of PSMs per protein per group 
                              # (Original, New, Total)
                              # to build the correlation matrix
) {

  library(tidyverse)
  library(reshape2)
  library(stringr)
  source('lib.R')
  
  #cat('Loading evidence file...\n')
  # load evidence file
  #ev <- parse.ev.adj(ev.in, type=ev.type)
  
  # assign protein IDs
  # ev$Protein_ID <- as.numeric(as.factor(ev$Proteins))
  
  cat('Validating raw files', paste(raw.files, collapse=', '), '\n')
  cat('Processing evidence file...\n')
  # # clean up raw file names and extract an experiment ID from it
  # # need the experiment ID (19A, 30B, etc.), to match with sample metadata from the
  # # experiment description excel sheet/.csv
  # ev$file <- clean.file.name(ev$`Raw file`)
  # ev$exp <- str_extract(ev$file, '[1-9][0-9][A-Z](\\d)?')
  # 
  # load excel description
  desc <- process.desc()
  desc.types <- dcast(desc, Exp~Type, value.var='Type', fun.aggregate=length)
  desc.types$Exp <- as.character(desc.types$Exp)
  
  if (!is.null(raw.files)) {
    # filter by raw files if defined and experiment IDs aren't
    if(length(raw.files) < 1) {
      stop('No raw file names given')
    }
    ev <- ev[ev$`Raw file` %in% raw.files,]
  } 
  
  # make sure at this point that we still have observations to work with
  if(nrow(ev) < 1) {
    stop('Incorrect or too severe filtering -- no data left')
  }
  
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
  
  if(remove.keratins) {
    # load keratin protein IDs
    keratins <- read_lines('Keratins.txt')
    # this is slow for many observations. might be faster to do a
    # str_extract to get the uniprot ID, and then to do a match() call
    ev <- ev[!grepl(paste(keratins, collapse='|'), ev$Proteins),]
  }
  
  # set zeroes to NA
  ev[,data.cols][ev[,data.cols]==0] <- NA
  
  # remove rows that have 50% or more of their columns as NAs
  num.nas <- apply(t(apply(ev[,data.cols], MARGIN=1, FUN=is.na)), MARGIN=1, FUN=sum)
  ev <- ev[ num.nas < ceiling(length(data.cols) * sparse.filter),]
  
  cat(nrow(ev), 'observations remaining after filtering\n')
  
  ## remove carrier and empty channels
  # these channel indices vary between experiments, 
  # so we're gonna have to loop thru these and check for each one
  exps <- str_extract(clean.file.name(unique(ev$`Raw file`)), '[1-9][0-9][A-Z]')
  
  for(i in unique(exps)) {
    #if(sum(grepl(i,ev$exp)) <= 0) next
    inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch']

    cat('Channel(s)',inds.remove,'from experiment',i,'are blanks or carrier channels. Removing...\n')
    
    # make relative to the ev data frame
    inds.remove <- inds.remove + data.cols[1] - 1
    # set to NA
    ev[,inds.remove] <- NA
  }
  
  ## normalize data
  
  # first normalize by column, by median
  # this assumes the same amount of protein per channel
  # (not by cell, but in aggregate - this erases any ionization or 
  # collision cell or any other quantitative biases)
  if(norm.cols) {
    ev[,data.cols] <- sweep(
      # operate on the data
      x=ev[,data.cols], 
      # operate over the columns (2)
      MARGIN=2, 
      # '/' specifies the division operator/function
      # can do '*', or any other function or overloaded operator)
      FUN='/',
      # divide by the median of each column (2). remove NAs
      STATS=apply(ev[,data.cols], MARGIN=2, FUN=median, na.rm=TRUE))
  }
  
  # now normalize across rows, to get the difference between the channels
  # AKA the difference between different cells/conditions
  if(norm.rows) {
    ev[,data.cols] <- sweep(
      # operate on the data
      x=ev[,data.cols],
      # operate over the rows (1)
      MARGIN=1, 
      # '/' specifies the division operator/function
      FUN='/',
      # divide by the median of each row (1). remove NAs
      STATS=as.vector(apply(ev[,data.cols], MARGIN=1, FUN=median, na.rm=TRUE)))
  }
  
  ## only use proteins with many peptides
  
  # count # of peptides per protein in our entire dataset
  prot.map <- as.data.frame(table(ev$Proteins))
  prot.map <- prot.map[order(prot.map$Freq, decreasing=TRUE),]
  
  # only use proteins with more than [pep.freq.thresh] peptides attached to it
  # pep.freq.thresh defaults to 10
  prots <- as.character(prot.map$Var1[prot.map$Freq >= psm.freq.thresh])
  
  if(length(prots) < 1) {
    stop(paste('No proteins with more PSMs than', psm.freq.thresh))
  }
  
  cat(length(prots), 'proteins selected to build correlation matrices\n')
  
  ## Correlate peptides within each protein
  #
  # Iterate over proteins and do pairwise correlation between 
  # quantitation vectors of peptides.
  # Do this for certain subsets of the protein data, 
  # in order to compare them after the fact
  
  cat('Building correlation vectors...\n')
  cors <- vector(mode='numeric')
  cors.new <- vector(mode='numeric')
  # cors.tot <- vector(mode='numeric')
  
  # significance threshold
  alpha <- 0.05
  
  #cor.size <- vector(mode='numeric')
  
  for(i in 1:length(prots)) {
    
    #prot <- 'sp|P06733|ENOA_HUMAN'
    #prot <- 'sp|P11142|HSP7C_HUMAN' # heat shock
    
    prot <- prots[i]
    cat('\r', match(prot, prots), '/', length(prots), '-', prot, '                           ')
    flush.console()
    
    # get PSMs for this protein
    
    # ORIGINAL 
    # -- originally good IDs (PEP < 0.05), ignoring result of RTLib method
    prot.data <- ev[ev$Proteins==prot & ev$PEP < alpha, data.cols]
    #prot.names <- ev[ev$Proteins==prot & ev$PEP < alpha,]$Sequence
    
    # NEW
    # -- originally bad IDs, upgraded to good ID (PEP > 0.05 & PEP.new < 0.05)
    # -- set of NEW and ORIGINAL should be disjoint, i.e., NEW ∩ ORIGINAL = ∅
    prot.data.new <- ev[ev$Proteins==prot & ev$PEP > alpha & ev$PEP.new < alpha, data.cols]
    #prot.names.new <- ev[ev$Proteins==prot & ev$PEP > alpha & ev$PEP.new < alpha,]$Sequence
    
    # TOTAL
    # -- good IDs after bayesian update, regardless of if it was upgraded beyond threshold or not
    # -- **NOT** the union of ORIGINAL and NEW (but does contain all of NEW), 
    # -- as some PSMs will be downgraded below PEP threshold via. the RTLib method
    # --
    # -- reminder: PEP.updated is PEP superceded by PEP.new, so that there are no NAs
    # -- that is, a PSM that was updated by the bayesian update will have its PEP replaced by PEP.new
    # prot.data.tot <- ev[ev$Proteins==prot & (ev$PEP.new < alpha | ev$PEP < alpha),data.cols]
    
    # init output
    prot.cor <- NULL
    prot.cor.new <- NULL
    # prot.cor.tot <- NULL
    
    # lets say... must have more than n observations to have a meaninful correlation matrix
    # not enough data can push extremely low/high correlations, especially when so many
    # of our columns are NA anyways
    #
    # the pairwise.complete.obs method of the cor() function can exacerbate 
    # this effect if there aren't enough observations
    #cor.size.thresh <- 5
    
    if(nrow(prot.data) < cor.size.thresh | nrow(prot.data.new) < cor.size.thresh) {
      next
    }
    
    # pairwise complete to account to missing/NA values
    # results in matrix that is not positive semi-definite, and rows that may contain NA
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    # set diagonal to NA
    diag(prot.cor) <- NA

    # do the same for prot.data.new, prot.data.tot
    prot.cor.new <- cor(t(prot.data.new), use='pairwise.complete.obs', method='pearson')
    diag(prot.cor.new) <- NA
      
    # prot.cor.tot <- cor(t(prot.data.tot), use='pairwise.complete.obs', method='pearson')
    # diag(prot.cor.tot) <- NA
    
    # concatenate to master vectors
    cors <- c(cors, as.vector(prot.cor))
    cors.new <- c(cors.new, as.vector(prot.cor.new))
    # cors.tot <- c(cors.tot, as.vector(prot.cor.tot))

    # cor.size <- c(cor.size, length(prot.cor))
  }
  
  #plot(cumsum(cor.size), xlab='Proteins', ylab='Correlation Vector Size')
  #abline(v=30, col='red')
  
  ## Compute null distribution for peptide/protein correlations
  
  cat('\nBuilding null correlation vector...\n')
  cors.null <- vector(mode='numeric')
  # set seed for consistency
  set.seed(1)
  # take n at a time
  # n = median # of observations for this set of proteins
  # this should prevent too many/too little observations from affecting
  # the null distribution vs. the correct ones
  # n <- median(prot.map$Freq[prot.map$Var1 %in% prots])
  
  for(i in 1:length(prots)) {
    cat('\r', i, '/', length(prots), '-', 'dummy', i, '                           ')
    flush.console()
    
    ## pairwise correlation between peptides from random proteins
    if(i < 1) next
    
    # randomly sample peptides from this list, and then compute the
    # pairwise correlation between them
    #
    # ideally this will be uniform - but will probably never be as
    # that would assume that all proteins are regulated independently
    # there are definitely proteins whose levels will correlate closely with others
    # expect this to be only slightly biased in the positive direction
    
    n <- prot.map[i,'Freq']
    
    prot.data <- ev[sample(1:nrow(ev), size=n), data.cols]
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    
    # set diagonal to NA
    diag(prot.cor) <- NA
    # concatenate to cors.null vector
    cors.null <- c(cors.null, as.vector(prot.cor))
    
    #plot(density(as.vector(prot.cor), na.rm=T))
  }
  
  # consolidate all groups into one dataframe
  cat('\nConsolidating correlation vectors...\n')
  cors.all <- data.frame(
    # data=as.numeric(c(cors, cors.new, cors.tot, cors.null)),
    data=as.numeric(c(cors, cors.new, cors.null)),
    Type=as.factor(c(rep('Original', length(cors)),
                     rep('New', length(cors.new)),
                     # rep('Total', length(cors.tot)),
                     rep('Null', length(cors.null))))
  )
  # remove NAs
  cors.all <- cors.all[!is.na(cors.all$data),]
  
  # reorder factors
  #cors.all$Type <- fct_relevel(cors.all$Type, 'Original', 'New', 'Null', 'Total')
  
  # plot(density(cors.new, na.rm=T), col='blue')
  # lines(density(cors, na.rm=T))
  # lines(density(cors.null, na.rm=T), col='red')
  
  
  return(cors.all)
}