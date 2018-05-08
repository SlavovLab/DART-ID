validate.lib <- function(
  ev, 
  ev.type='STAN',
  #exps=NULL,                 # subset of experiments to use
  raw.files=NULL,             # subset of raw files to use
  sparse.filter=0.5,          # filter out sparse observations. 
                              # this is the fraction of columns in 
                              # the TMT data that is allowed to be NA
  unique.only=TRUE,           # only take unique peptides
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
  
  cat('Validating raw files', paste(raw.files, collapse=', '), '\n')
  cat('Processing evidence file...\n')

  # load excel description
  desc <- process.desc()
  desc.types <- dcast(desc, Exp~Type, value.var='Type', fun.aggregate=length)
  desc.types$Exp <- as.character(desc.types$Exp)
  
  if (!is.null(raw.files)) {
    # filter by raw files if defined and experiment IDs aren't
    if(length(raw.files) < 1) {
      #stop('No raw file names given')
    }
    ev <- ev[ev$`Raw file` %in% raw.files,]
  } 
  
  # find out the indices of columns with TMT reporter ion data
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
  
  # remove REV and CON proteins
  ev <- ev[!grepl('REV*',ev$Proteins),]
  ev <- ev[!grepl('CON*',ev$Proteins),]
  
  # load keratin protein IDs
  keratins <- read_lines('Keratins.txt')
  # this is slow for many observations. might be faster to do a
  # str_extract to get the uniprot ID, and then to do a match() call
  ev <- ev[!grepl(paste(keratins, collapse='|'), ev$Proteins),]
  
  # set zeroes to NA
  ev[,data.cols][ev[,data.cols] == 0] <- NA
  
  # remove rows with more than (sparse.filter as fraction) NAs
  ev <- ev[apply(!apply(ev[,data.cols], 1, is.na), 2, sum) >= length(data.cols) * sparse.filter,]
  
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
  ev[,data.cols] <- t(t(ev[,data.cols]) / apply(ev[,data.cols], 2, median, na.rm=T))
  
  # now normalize across rows, to get the difference between the channels
  # AKA the difference between different cells/conditions
  ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], 1, median, na.rm=T)
  #ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], 1, mean, na.rm=T)
  
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
  
  cat(length(prots), 'proteins with over', psm.freq.thresh, 'PSMs selected to build correlation matrices\n')
  
  ## Correlate peptides within each protein
  #
  # Iterate over proteins and do pairwise correlation between 
  # quantitation vectors of peptides.
  # Do this for certain subsets of the protein data, 
  # in order to compare them after the fact
  
  cat('Building correlation vectors...\n')
  cors <- vector(mode='numeric')
  cors.new <- vector(mode='numeric')
  
  # significance threshold
  alpha <- 0.05
  
  for(i in 1:length(prots)) {
    
    #prot <- 'sp|P06733|ENOA_HUMAN'
    #prot <- 'sp|P11142|HSP7C_HUMAN' # heat shock
    
    prot <- prots[i]
    cat('\r', match(prot, prots), '/', length(prots), '-', prot, '                           ')
    flush.console()
    
    # get PSMs for this protein
    
    # ORIGINAL 
    # -- originally good IDs (PEP < 0.05), ignoring result of RTLib method
    prot.data <- data.matrix(subset(ev, Proteins==prot & PEP < alpha, 
                                    select=data.cols))
    prot.names <- ev[ev$Proteins==prot & ev$PEP < alpha,] %>%
      pull('Sequence')
    
    # NEW
    # -- originally bad IDs, upgraded to good ID (PEP > 0.05 & PEP.new < 0.05)
    # -- set of NEW and ORIGINAL should be disjoint, i.e., NEW ∩ ORIGINAL = ∅
    prot.data.new <- data.matrix(subset(
      ev, Proteins==prot & PEP > alpha & PEP.new < alpha, select=data.cols))
    prot.names.new <- subset(ev, Proteins==prot & PEP > alpha & PEP.new < alpha) %>%
      pull('Sequence')
    
    # lets say... must have more than n observations to 
    # have a meaninful correlation matrix
    # not enough data can push extremely low/high correlations, 
    # especially when so many of our columns are NA anyways
    #
    # the pairwise.complete.obs method of the cor() function can exacerbate 
    # this effect if there aren't enough observations
    #cor.size.thresh <- 5
    
    if(nrow(prot.data) < cor.size.thresh | 
       nrow(prot.data.new) < 1) {
      next
    }
    if(similarity=='correlation') {
      prot.cor <- cor.sim(prot.data)
      prot.cor.new <- cor.sim(rbind(prot.data.new, prot.data))
    } else if (similarity == 'cosine') {
      prot.cor <- cos.sim(prot.data)
      prot.cor.new <- cos.sim(prot.data)
    }
    
    cor.mat <- prot.cor.new
    #cor.mat <- prot.cor
    diag(cor.mat) <- 1
    #cor.mat[upper.tri(cor.mat)] <- NA
    rownames(cor.mat) <- c(prot.names.new, prot.names)
    colnames(cor.mat) <- c(prot.names.new, prot.names)
    # rownames(cor.mat) <- prot.names
    # colnames(cor.mat) <- prot.names
    
    # ggplot(melt(cor.mat, na.rm=T), aes(x=Var1, y=Var2, fill=value)) +
    #   geom_tile(color='white') +
    #   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    #                        midpoint = 0, limit = c(-1,1), space = "Lab", 
    #                        name="Pearson\nCorrelation") +
    #   theme_minimal()+ 
    #   theme(axis.text.y = element_text(size=8),
    #         axis.text.x = element_text(angle = 45, vjust = 1, 
    #                                    size = 8, hjust = 1)
    #         #axis.text.x = element_blank()
    #         ) +
    #   labs(x=NULL, y=NULL)+
    #   coord_fixed()
    
    # concatenate to master vectors
    cors <- c(cors, as.vector(prot.cor))
    # only take comparisons between new data itself and new data vs. old data,
    # don't copy the original correlation matrix in
    cors.new <- c(cors.new, c(
      as.vector(prot.cor.new[1:nrow(prot.data.new),-(1:nrow(prot.data.new))]),
      as.vector(prot.cor.new[,1:nrow(prot.data.new)])))

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
    
    if(similarity=='correlation') {
      prot.cor <- cor.sim(prot.data)
    } else if (similarity=='cosine') {
      prot.cor <- cos.sim(prot.data)
    }
    
    # set diagonal to NA
    diag(prot.cor) <- NA
    # concatenate to cors.null vector
    cors.null <- c(cors.null, as.vector(prot.cor))
    
    #plot(density(as.vector(prot.cor), na.rm=T))
  }
  
  # remove NAs
  cors <- cors[!is.na(cors)]
  cors.new <- cors.new[!is.na(cors.new)]
  
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
  
  #cor <- cor(t(data), method='pearson')
  cor <- cor(t(data), method='spearman')
  # set diagonal to NA
  diag(cor) <- NA
  
  return(cor)
}
