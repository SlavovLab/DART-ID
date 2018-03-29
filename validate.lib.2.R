validate.lib <- function(
  ev, 
  sparse.filter=0.5,          # filter out sparse observations. 
                              # this is the fraction of columns in 
                              # the TMT data that is allowed to be NA
  unique.only=FALSE,          # only take unique peptides
  psm.freq.thresh=10,         # minimum number of PSMs per protein
  alpha=0.05,                 # significance threshold that separates
                              # Original and New sets (PEP < alpha)
  cor.size.thresh=5,          # minimum number of PSMs per protein per group 
                              # (Original, New, Total)
                              # to build the correlation matrix
  similarity='correlation'    # 'correlation' for pearson correlation
                              # 'cosine' for cosine similarity
) {

  library(tidyverse)
  library(reshape2)
  library(stringr)
  # for imputation
  library(VIM)
  source('lib.R')
  
  #cat('Loading evidence file...\n')
  # load evidence file
  #ev <- parse.ev.adj(ev.in, type=ev.type)
  
  # assign protein IDs
  # ev$Protein_ID <- as.numeric(as.factor(ev$Proteins))
  
  #cat('Validating raw files', paste(raw.files, collapse=', '), '\n')
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
  ev <- ev[!grepl('REV*',ev$Proteins),]
  
  # remove CON proteins
  # common contaminants, like Keratin, Albumin, etc.
  ev <- ev[!grepl('CON*',ev$Proteins),]
  
  # load keratin protein IDs
  keratins <- read_lines('Keratins.txt')
  # this is slow for many observations. might be faster to do a
  # str_extract to get the uniprot ID, and then to do a match() call
  ev <- ev[!grepl(paste(keratins, collapse='|'), ev$Proteins),]

  # set zeroes to NA
  ev[,data.cols][ev[,data.cols]==0] <- NA
  
  ## remove carrier and empty channels
  # these channel indices vary between experiments, 
  # so we're gonna have to loop thru these and check for each one
  #exps <- str_extract(clean.file.name(unique(ev$`Raw file`)), '[1-9][0-9][A-Z]')
  
  # pull experiment number from raw file name
  ev$exp <- str_extract(clean.file.name(ev$`Raw file`), '[1-9][0-9][A-Z]')
  # remove rows that don't belong to a recognized experiment
  ev <- ev[!is.na(ev$exp),]
  
  for(i in unique(ev$exp)) {
    #if(sum(grepl(i,ev$exp)) <= 0) next
    inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch']

    cat('Channel(s)',inds.remove,'from experiment',i,'are blanks or carrier channels. Removing...\n')
    
    # make relative to the ev data frame
    inds.remove <- inds.remove + data.cols[1] - 1
    # set to NA
    ev[ev$exp==i,inds.remove] <- NA
  }
  
  cat(nrow(ev), 'observations remaining after filtering\n')
  
  # remove rows with more than (sparse.filter as fraction) NAs
  ev <- ev[apply(!apply(ev[,data.cols], 1, is.na), 2, sum) >= length(data.cols) * sparse.filter,]
  
  ## normalize data
  
  # first normalize by column, by median
  # this assumes the same amount of protein per channel
  # (not by cell, but in aggregate - this erases any ionization or 
  # collision cell or any other quantitative biases)
  ev[,data.cols] <- t(t(ev[,data.cols]) / apply(ev[,data.cols], 2, median, na.rm=T))
  
  # now normalize across rows, to get the difference between the channels
  # AKA the difference between different cells/conditions
  #ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], 1, median, na.rm=T)
  ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], 1, mean, na.rm=T)
  
  ## only use proteins with many peptides
  
  # count # of peptides per protein in our entire dataset
  prot.map <- as.data.frame(table(ev$Proteins))
  prot.map <- prot.map[order(prot.map$Freq, decreasing=TRUE),]
  
  # only use proteins with more than [pep.freq.thresh] peptides attached to it
  # pep.freq.thresh defaults to 10
  prots <- as.character(prot.map$Var1[prot.map$Freq >= psm.freq.thresh])
  
  cat(length(prots), 'proteins selected to build correlation matrices\n')
  
  ## Correlate peptides within each protein
  # Iterate over proteins and do pairwise correlation between 
  # quantitation vectors of peptides.
  # Do this for certain subsets of the protein data, 
  # in order to compare them after the fact
  
  cat('Building correlation vectors...\n')
  cors <- vector(mode='numeric')
  cors.new <- vector(mode='numeric')
  
  # significance threshold
  alpha <- 0.05
  raw.files <- unique(ev$`Raw file`)
  
  for(i in 1:length(prots)) {
    
    prot <- prots[i]
    cat('\r', match(prot, prots), '/', length(prots), '-', prot, '                           ')
    flush.console()
    
    # get PSMs for this protein
    ev.p <- subset(ev, Proteins==prot & PEP < alpha)
    ev.p.new <- subset(ev, Proteins==prot & PEP > alpha & PEP.new < alpha)
    
    if(nrow(ev.p) < 1 & nrow(ev.p.new) < 1) {
      next
    }
    
    #prot.data <- ev.p[,data.cols]
    prot.data <- matrix(nrow=(nrow(ev.p) + nrow(ev.p.new)), 
                        ncol=length(data.cols)*length(raw.files))
    
    row.counter <- 1
    row.counter.new <- nrow(prot.data) - nrow(ev.p.new) - 1
    if(row.counter.new < 1) { row.counter.new <- 1 }
    
    for(j in 1:length(raw.files)) {
      raw.file <- raw.files[j]
      
      ev.pr <- data.matrix(subset(ev.p, `Raw file`==raw.file, select=data.cols))
      if(nrow(ev.pr) > 0) {
        prot.data[row.counter:(row.counter+nrow(ev.pr) - 1),
                  (((j - 1) * length(data.cols)) + 1):(j * length(data.cols))] <- ev.pr
        row.counter <- row.counter + nrow(ev.pr)
      }
      
      ev.pr.new <- data.matrix(subset(ev.p.new, `Raw file`==raw.file, select=data.cols))
      if(nrow(ev.pr.new) > 0) { 
        prot.data[row.counter.new:(row.counter.new + nrow(ev.pr.new) - 1),
                  (((j - 1) * length(data.cols)) + 1):(j* length(data.cols))] <- ev.pr.new
        row.counter.new <- row.counter.new + nrow(ev.pr.new)
      }
    }
    
    # remove NA columns
    prot.data <- prot.data[,apply(prot.data, 2, sum, na.rm=T)!=0]
    
    if(length(prot.data) < 1) next
    
    # run spearman correlation on prot data matrix
    # this might take a while...
    cor.mat <- cor(t(prot.data), use='pairwise.complete.obs', method='spearman')
    # eliminate diagonal
    diag(cor.mat) <- NA
    
    # squash into 1-d vectors
    prot.cor <- as.numeric(cor.mat[1:nrow(ev.p),1:nrow(ev.p)])
    prot.cor <- prot.cor[!is.na(prot.cor)]
    if(nrow(ev.p.new) > 0) {
      prot.cor.new <- as.numeric(cor.mat[(nrow(ev.p)+1):nrow(cor.mat),])
      prot.cor.new <- prot.cor.new[!is.na(prot.cor.new)]
    }
    
    # concatenate to master vectors
    cors <- c(cors, prot.cor)
    cors.new <- c(cors.new, prot.cor.new)
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
    
    # randomly sample peptides from this list, and then compute the
    # pairwise correlation between them
    #
    # ideally this will be uniform - but will probably never be as
    # that would assume that all proteins are regulated independently
    # there are definitely proteins whose levels will correlate closely with others
    # expect this to be only slightly biased in the positive direction
    
    n <- prot.map[i,'Freq']
    ev.r <- ev[sample(1:nrow(ev), size=n),]
    prot.data <- matrix(nrow=n, ncol=length(data.cols) * length(raw.files))
    
    row.counter <- 1
    for(j in 1:length(raw.files)) {
      raw.file <- raw.files[j]
      
      ev.rr <- data.matrix(subset(ev.r, `Raw file`==raw.file, select=data.cols))
      if(nrow(ev.rr) > 0) {
        prot.data[row.counter:(row.counter+nrow(ev.rr) - 1),
                  (((j - 1) * length(data.cols)) + 1):(j * length(data.cols))] <- ev.rr
        row.counter <- row.counter + nrow(ev.rr)
      }
    }
    
    # remove NA columns
    prot.data <- prot.data[,apply(prot.data, 2, sum, na.rm=T)!=0]
    
    # run spearman correlation on prot data matrix
    # this might take a while...
    cor.mat <- cor(t(prot.data), use='pairwise.complete.obs', method='spearman')
    
    # set diagonal to NA
    diag(cor.mat) <- NA
    
    # 1-dimensionalise correlation matrix
    null.cor <- as.vector(cor.mat)
    null.cor <- null.cor[!is.na(null.cor)]
    
    # concatenate to cors.null vector
    cors.null <- c(cors.null, null.cor)
    
    #plot(density(as.vector(prot.cor), na.rm=T))
  }
  
  # par(mfrow=c(1,2))
  # plot(x=seq(-1,1,by=0.01), y=1-ecdf(cors)(seq(-1,1,by=0.01)), type='l', col='black',
  #      xlab='Correlation', ylab='Survival (1-CDF)', main='Validation - Mean Normalization')
  # lines(x=seq(-1,1,by=0.01), y=1-ecdf(cors.new)(seq(-1,1,by=0.01)), type='l', col='blue')
  # lines(x=seq(-1,1,by=0.01), y=1-ecdf(cors.null)(seq(-1,1,by=0.01)), type='l', col='red')
  # 
  # plot(density(cors), col='black', xlab='Correlation', ylab='Density', main='')
  # lines(density(cors.new), col='blue')
  # lines(density(cors.null), col='red')
  
  save(cors, cors.new, cors.null, file='cors.validation.2.mean.RData')
  
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

cos.sim <- function(data) {
  ## cosine angle between vectors of RI intensity, instead of correlation
  # remove NA cols
  data <- data[,!apply(!apply(data, 1, is.na), 1, sum) == 0]
  # remove NA rows
  data <- data[!apply(!apply(data, 1, is.na), 2, sum) == 0,]
  # impute missing data
  data <- kNN(data, k=5, imp_var=F)
  # normalize by vector length
  data <- data.matrix(data / apply(data, 1, FUN=function(x) {
    sqrt(sum(x^2))
  }))
  # compute cosine similarity
  # this is the same as the normalized dot product,
  # but since normalization was already done, this can simply be done
  # by doing A * At
  cor <- data %*% t(data)
  
  return(cor)
}

cor.sim <- function(data) {
  # remove NA cols
  data <- data[,!apply(!apply(data, 1, is.na), 1, sum) == 0]
  # remove NA rows
  data <- data[!apply(!apply(data, 1, is.na), 2, sum) == 0,]
  # impute missing data
  data <- kNN(data, k=5, imp_var=F)
  
  cor <- cor(t(data), method='pearson')
  # set diagonal to NA
  diag(cor) <- NA
  
  return(cor)
}
