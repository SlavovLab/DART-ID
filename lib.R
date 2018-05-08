source('validate.lib.R')

load.ev <- function(ev, par.file='dat/params.Fit2.RData', include.REV=FALSE, include.CON=FALSE) {
  load(par.file)
  # remove abnormal LC experiments
  # load experiments from correlation testing in similar.lc.R
  exps.lc <- unlist(read_csv('dat/exps.corr.txt')[,2])
  names(exps.lc) <- NULL
  
  ev <- ev %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>% # Only use Elite experiments
    filter(`Raw file` %in% exps.lc) # Remove abnormal LC experiments
  
  ## Filter of PEP < .05
  ev.f <- ev %>% filter(PEP < 0.05) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) # Only use Elite experiments

  # Remove Reverse matches
  if(!include.REV) {
    ev.f <- ev.f %>% filter(!grepl('REV*', `Leading razor protein`)) 
  }  
  # Remove Contaminants
  if(!include.CON) {
    ev.f <- ev.f %>% filter(!grepl('CON*',`Leading razor protein`))
  }
  
  ev.f <- ev.f %>% 
    filter(PEP < 0.05) %>%
    mutate(Protein=`Leading razor protein`) %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP", "Protein") %>%
    mutate(exp_id=`Raw file`) %>%  # new column - exp_id = numeric version of experiment file
    mutate_at("exp_id", funs(as.numeric(as.factor(.)))) %>%
    mutate(`Stan ID`=`Peptide ID`) %>%
    mutate_at("Stan ID", funs(as.numeric(as.factor(.))))
  ev.f <<- ev.f
  
  experiment_factors <<- as.factor(ev.f$`Raw file`)
  num_exps <<- length(unique(ev.f[["exp_id"]]))
  
  ## "true peptide id" matches peptide id in evidence file
  ## "peptide id" is index in 1:num_peptides for stan
  raw_peptide_id <<- ev.f[["Peptide ID"]]
  pep.id.list <<- unique(raw_peptide_id)
  num_peptides <<- length(unique(ev.f$`Stan ID`))
  
  # parse linear regression params from STAN output
  beta0 <<- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
  beta1 <<- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
  beta2 <<- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
  
  # sigma params from Fit2
  split.point <<- pars[sprintf('split_point[%i]', seq(1, num_exps))]
  sigma.slope <<- pars[sprintf('sigma_slope[%i]', seq(1, num_exps))]
  sigma.intercept <<- pars[sprintf('sigma_intercept[%i]', seq(1, num_exps))]
  sigma.slope.global <<- pars['sigma_slope_global']
  
  # sigma params from Fit1
  sigmas <<- pars[grep('sigma\\[', names(pars))]
  
  mus <<- pars[grep('mu\\[', names(pars))]
}

clean.file.name <- function(name) {
  name <- gsub('#', '', name)
  name <- gsub('.*_NC_', '', name)
  name <- gsub('.raw', '', name)
  name <- gsub('set', '', name)
  name <- gsub('[-|+|=]', '_', name)
  return(name)
}

theme_bert <- function() {
  theme_bw(base_size=16, base_family="Helvetica") %+replace% 
    theme(
      #panel.background  = element_rect(fill=NULL, color='black', size=0.5),
      panel.background=element_blank(),
      panel.border=element_blank(),
      axis.line.x=element_line(color='black', size=0.25),
      axis.line.y=element_line(color='black', size=0.25),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      #axis.line=element_line(size=0.5, color='#888888'),
      #axis.line=element_blank(),
      axis.ticks=element_line(color='black', size=0.25),
      #axis.ticks=element_blank(),
      panel.grid=element_blank()
    )
}

fold.change.comp <- function(exps, begin=1e-5, end=1, num.steps=100, log=T, add.dummy=T) {
  library(pracma)
  # only use the same raw files between all evidence files
  common.exps <- Reduce(intersect, lapply(exps, function(exp) { exp$`Raw file` }))
  exps <- lapply(exps, function(exp) { exp[exp$`Raw file` %in% common.exps,]})
  
  # equally spaced steps in log space
  if(log) {
    x <- logseq(begin, end, n=num.steps)
  } else {
    x <- seq(begin, end, length.out=num.steps)
  }
  
  # frame to hold the results
  df <- data.frame()
  counter <- 0
  for(i in x) {
    ratios <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.new < i, na.rm=TRUE) / 
         sum(exp$PEP < i & !is.na(exp$PEP.new)))
    }))
    
    # add dummy maxquant experiment
    if(add.dummy) {
      ratios <- c(ratios, 1)
      exp.names <- c(names(exps), 'MaxQuant')
    } else {
      exp.names <- names(exps)
    }
    
    df <- rbind(df, data.frame(
      x=as.numeric(i),
      PEP=as.numeric(ratios),
      Method=as.character(exp.names)
    ))
  }
  return(df)
}

# fancy scientific scales
# from: https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  #l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  #l <- gsub("e", "%*%10^", l)
  # make sure +0 just turns into 0
  l <- gsub("\\+00", "00", l)
  # return this as an expression
  return(parse(text=l))
}

# from: https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
emptyPlot <- function() {
  ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
}

bHeatmap <- c('#710e0d', '#8b42df', '#4c94dc', '#40df91', '#fef245')

#' Process Experiment Description Excel Sheet
process.desc <- function() {
  desc <- read.csv('dat/SingleCellExperiments_Description.csv')
  # only take the portion we need
  desc <- desc[1:308,]
  # remove rows that don't have experiment info
  desc <- desc[grepl('^\\d{2,3}[A-Z]{1}', desc$`Exp...`),]
  desc <- desc[,c(1, 3:12)]
  names(desc)[1] <- 'Exp'
  rownames(desc) <- NULL
  desc <- melt(desc, id.vars=c('Exp'), variable.name='Channel', value.name='Sample',
               factorsAsStrings=F)
  desc$Exps <- as.character(desc$Exp)
  # channel ID
  desc$ch <- as.numeric(desc$Channel)
  
  desc$Sample <- as.character(desc$Sample)
  desc$Sample[desc$Sample %in% c('', '0', 'Empty', 'N/A', 'PBS')] <- NA
  
  ## Quantity - how many cells were in each channel
  desc$Quantity <- str_extract(desc$Sample, '\\d+(e\\d)?(\\.\\dk?)?')
  # convert scientific notation
  desc$Quantity[grep('\\de\\d', desc$Quantity)] <- 
    as.numeric(desc$Quantity[grep('\\de\\d', desc$Quantity)])
  # convert "k" notation
  desc$Quantity[grep('\\d\\.\\dk', desc$Quantity)] <- 
    as.numeric(str_extract(desc$Quantity[grep('\\d+\\.\\dk', desc$Quantity)], '\\d+\\.\\d')) * 1000
  
  ## Type - what kind of cell it was
  # J = Jurkat
  # H = Hek293
  # U = U937
  # ES|EB = Mouse
  desc$Type <- str_extract(desc$Sample, '[JHUjhu]{1}|([Ee]{1}[SsBb]{1})')
  # Mixed type - more than one type of cell
  # most likely a carrier channel
  desc$Type[grepl('\\&|and', desc$Sample)] <- 'Mixed'
  
  # Uppercase = intact, single cell
  # Lowercase = diluted, lysed mixture
  desc$Diluted <- F
  desc$Diluted[grepl('[jhu]|es|eb', desc$Type)] <- T
  desc$Diluted[is.na(desc$Type)] <- NA
  
  # Re-normalize cell type for easier searching
  desc$Type <- toupper(desc$Type)
  
  return(desc)
}

fold.change.comp <- function(exps) {
  library(pracma)
  # only use the same raw files between all evidence files
  #common.exps <- Reduce(intersect, lapply(exps, function(exp) { exp$`Raw file` }))
  #exps <- lapply(exps, function(exp) { exp[exp$`Raw file` %in% common.exps,]})
  
  num.steps <- 100
  # equally spaced steps in log space
  #x <- logseq(1e-10, 1, n=num.steps)
  x <- logseq(1e-5, 1, n=num.steps)
  # frame to hold the results
  df <- data.frame()
  counter <- 0
  for(i in x) {
    ratios <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.new < i, na.rm=TRUE) / 
         sum(exp$PEP < i & !is.na(exp$PEP.new)))
    }))
    df <- rbind(df, data.frame(
      x=as.numeric(i),
      PEP=as.numeric(ratios),
      Method=as.character(names(exps))
    ))
  }
  return(df)
}

remove.channels <- function(ev) {
  # load excel description
  desc <- process.desc()
  desc.types <- dcast(desc, Exp~Type, value.var='Type', fun.aggregate=length)
  desc.types$Exp <- as.character(desc.types$Exp)
  
  ## remove carrier and empty channels
  # these channel indices vary between experiments, 
  # so we're gonna have to loop thru these and check for each one
  exps <- str_extract(clean.file.name(unique(ev$`Raw file`)), '[1-9][0-9][A-Z]')
  exps <- unique(exps[!is.na(exps)])
  
  for(i in exps) {
    #if(sum(grepl(i,ev$exp)) <= 0) next
    inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch']
    
    cat('Channel(s)',inds.remove,'from experiment',i,'are blanks or carrier channels. Removing...\n')
    
    # make relative to the ev data frame
    inds.remove <- inds.remove + data.cols[1] - 1
    # set to NA
    ev[,inds.remove] <- NA
  }
  
  return(ev)
}

get.exp.ids <- function(raw.files) {
  exps <- str_extract(clean.file.name(raw.files), '[1-9][0-9][A-Z]')
  exps <- unique(exps[!is.na(exps)])
  return(exps)
}

filter.exps <- function(ev, pep.thresh=1e-2) {
  # ev <- read_tsv('dat/evidence_elite.txt')
  
  #removed.exps <- c()
  
  # confidence threshold - harsh to increase confidence
  # and decrease computation time
  # pep.thresh <- 1e-2
  
  # make sure we have enough observations in all experiments
  # lets say... need at least 100 to get a good idea of experiment quality
  #below.thresh <- ev %>% 
  #  group_by(`Raw file`) %>% 
  #  summarise(n=sum(PEP < pep.thresh)) %>%
  #  arrange(desc(n)) %>%
  #  filter(n < 100) %>%
  #  pull(`Raw file`)
  
  #ev <- ev %>% filter(!`Raw file` %in% below.thresh)
  
  raw.files <- unique(ev$`Raw file`)
  cor.mat <- matrix(nrow=length(raw.files), ncol=length(raw.files))
  
  # correlate high-confidence PSMs between experiments
  # only need the upper triangular of the correlation matrix
  for(i in 1:length(raw.files)) {
    ev.a <- ev %>% 
      filter(PEP < pep.thresh) %>%
      filter(`Raw file`==raw.files[i]) %>%
      select(c('Mod. peptide ID', 'Retention time'))
    
    
    for(j in i:length(raw.files)) {
      if(i == j) next
      
      cat('\r', i, '-', j, '       ')  
      flush.console()
      
      ev.b <- ev %>%
        filter(PEP < pep.thresh) %>%
        filter(`Raw file`==raw.files[j]) %>%
        filter(`Mod. peptide ID` %in% ev.a$`Mod. peptide ID`) %>%
        select(c('Mod. peptide ID', 'Retention time')) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`))
      
      # require at least 5 common points to do this analysis
      if(nrow(ev.b) < 5) {
        cor.mat[i,j] <- NA
        next
      }
      
      b <- ev.b$`Retention time`
      a <- ev.a %>%
        filter(`Mod. peptide ID` %in% ev.b$`Mod. peptide ID`) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`)) %>%
        pull(`Retention time`)
      
      cor.mat[i,j] <- cor(a, b)
    }
  }
  
  return(cor.mat)
}

plot.cor.mat <- function(cor.mat, show.text=F) {
  library(ggplot2)
  library(reshape2)
  
  # only take upper triangle. don't want to destroy the matrix structure,
  # so set the lower triangle to a dummy value so it can be removed later
  cor.mat[lower.tri(cor.mat)] <- 100
  
  # set diagnoal to 1
  diag(cor.mat) <- 1
  
  df <- melt(cor.mat)
  df <- df %>%
    filter(!value == 100)
  
  p <- ggplot(df, aes(y=Var1, x=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value='grey50',
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(axis.text.y = element_text(size=8),
          axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 8, hjust = 1)
          #axis.text.x = element_blank()
    ) +
    labs(x=NULL, y=NULL) +
    coord_fixed()
  
  if(show.text) {
    p <- p + geom_text(aes(label=format(value, digits=2)))
  }
  
  p
}