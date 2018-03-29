library(tidyverse)
library(VIM)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.txt')

exps <- c('set#37A', 'set#37B', 'set#37C', 'set#37D',
          'set#38A', 'set#38B', 'set#38C', 'set#38D', 'set#39A', 'set#39B',
          'set#40A', 'set#40B', 'set#40C', 'set#40D',
          'set#32B', 'set#32C', 'set#32D',
          paste0('set#30', LETTERS[1:10]), 'set#29B',
          'set#29A', 'set#29C',
          paste0('set24', c('A','B','C')),
          paste0('set25', c('A','B','C')),
          paste0('set#26', c('A','B','C','D','E')),
          paste0('set#27', c('A','B','C','D')),
          paste0('set#28', c('A','B','C','D')),
          'set19B', '160712A_NC_19A',
          'set#30', 'set#29', '160822A_NC_set#27', '160913A_NC_set#32',
          '160712A_NC_19A', '160720A_NC_19A', '160727A_NC_19A', 
          '160801A_NC_set19A', '160808A_NC_set19A',
          '160818A_NC_set#19A', '160822A_NC_set#19A', '160831A_NC_set#19A',
          '160927A_NC_set#19A')
all.raw.files <- unique(ev$`Raw file`)

# find out the indices of columns with TMT reporter ion data
# should be 10 columns, but the code doesn't rely on this number
data.cols <- grep('intensity', colnames(ev))

# remove peptides without a protein
# this can happen, especially when searching with a restricted FASTA file
ev <- ev[!is.na(ev$Proteins),]

# remove REV and CON proteins
ev <- ev[!grepl('REV*',ev$Proteins),]
ev <- ev[!grepl('CON*',ev$Proteins),]

# load keratin protein IDs
keratins <- read_lines('Keratins.txt')
# this is slow for many observations. might be faster to do a
# str_extract to get the uniprot ID, and then to do a match() call
ev <- ev[!grepl(paste(keratins, collapse='|'), ev$Proteins),]

# load excel description
desc <- process.desc()
desc.types <- dcast(desc, Exp~Type, value.var='Type', fun.aggregate=length)
desc.types$Exp <- as.character(desc.types$Exp)

for(i in exps) {
  cat(i, '\n')
  
  raw.files <- all.raw.files[grep(i, all.raw.files)]
  
  if (is.null(raw.files)) { next } 
  
  exp <- first(get.exp.ids(raw.files))
  types <- desc[desc$Exp==exp & !(desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'Type']
  
  ev.a <- ev[ev$`Raw file` %in% raw.files,] %>% filter(PEP < 0.05)
  ev.a[,desc[desc$Exp==exp & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch'] + data.cols[1] - 1] <- NA
  dat.a <- data.matrix(ev.a[,data.cols])
  a <- go.pca(dat.a)
  
  ev.b <- ev[ev$`Raw file` %in% raw.files,] %>% filter(PEP > 0.05) %>% filter(PEP.new < 0.05)
  ev.b[,desc[desc$Exp==exp & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch'] + data.cols[1] - 1] <- NA
  dat.b <- data.matrix(ev.b[,data.cols])
  b <- go.pca(dat.b)
  
  df <- data.frame(
    x=as.numeric(c(a[1,],b[1,])),
    y=as.numeric(c(a[2,],b[2,])),
    type=as.character(rep(types,2)),
    set=as.character(rep(c('Original', 'New'),each=ncol(a)))
  )
  
  p <- 
  ggplot(df, aes(x=x, y=y, color=type, alpha=set)) +
    geom_hline(yintercept=0, color='black') +
    geom_vline(xintercept=0, color='black') +
    geom_point(size=6) +
    geom_text(aes(label=rep(seq(1,ncol(a)), 2)), color='white') +
    scale_alpha_manual(values=c(1,0.6), guide=F) +
    scale_color_manual(values=c('red', 'blue')) +
    labs(x='PC1', y='PC2') +
    theme_bert()
  
  ggsave(filename=paste('pca_plots/',i,'.pdf', sep=''), plot=p, device=pdf)
}


go.pca <- function(data) {
  # set 0s to NAs
  data[data == 0] <- NA
  # remove NA cols
  data <- data[,!apply(!apply(data, 1, is.na), 1, sum) == 0]
  # remove rows with more than 50% NA
  data <- data[apply(!apply(data, 1, is.na), 2, sum) > 5,]
  # normalize columns by median
  #data <- sweep(x=data, MARGIN=2, FUN='/',
  #              STATS=apply(data, MARGIN=2, FUN=median, na.rm=TRUE))
  # rescale rows around 0
  #data <- t(apply(data, 1, scale))
  # impute missing data
  data <- data.matrix(kNN(data, k=10, imp_var=F))
  if(sum(is.na(data)) > 0) {
    data <- data.matrix(kNN(data, k=10, imp_var=F))
  }
  
  # create similarity matrix
  sim <- t(data) %*% data
  # singular value decomposition
  pars <- svd(sim)
  
  #par(mfrow=c(1,1))
  #plot(pars$u[1,], pars$u[2,], xlab='PC1', ylab='PC2')
  #abline(h=0)
  #abline(v=0)
  
  return(pars$u)
}
