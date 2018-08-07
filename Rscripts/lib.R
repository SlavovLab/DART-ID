#source('validate.lib.R')

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
  theme_bw(base_size=10, base_family="Helvetica") %+replace% 
  #theme_bw(base_size=16, base_family="Helvetica") %+replace% 
    theme(
      #panel.background  = element_rect(fill=NULL, color='black', size=0.5),
      #panel.background=element_blank(),
      #panel.border=element_blank(),
      #axis.line.x=element_line(color='black', size=0.25),
      #axis.line.y=element_line(color='black', size=0.25),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      #axis.line=element_line(size=0.5, color='#888888'),
      axis.line=element_blank(),
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
vc <- c("#004358", "#1F8A70", "#BEDB39", "#FFE11A", "#FD7400")
av <- c("#20B2CF", "#85DB86", "#F29F05", "#F25C05", "#D92525")

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

fold.change.comp <- function(exps, range=c(1e-3, 1)) {
  library(pracma)
  # only use the same raw files between all evidence files
  #common.exps <- Reduce(intersect, lapply(exps, function(exp) { exp$`Raw file` }))
  #exps <- lapply(exps, function(exp) { exp[exp$`Raw file` %in% common.exps,]})
  
  num.steps <- 100
  # equally spaced steps in log space
  #x <- logseq(1e-10, 1, n=num.steps)
  x <- logseq(range[1], range[2], n=num.steps)
  # frame to hold the results
  df <- data.frame()
  counter <- 0
  for(i in x) {
    ratios <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.updated < i) / 
         sum(exp$PEP < i))
    }))
    identified_spectra <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP < i)) / nrow(exp)
    }))
    identified_spectra_new <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.updated < i)) / nrow(exp)
    }))
    df <- rbind(df, data.frame(
      x=as.numeric(i),
      PEP=as.numeric(ratios),
      ident=as.numeric(identified_spectra),
      ident_new=as.numeric(identified_spectra_new),
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

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

normalize_ri_data <- function(dmat) {
  dmat[dmat==0] <- NA
  
  ## normalize data
  
  # first normalize by column, by median
  dmat <- t(t(dmat) / apply(dmat, 2, median, na.rm=T))
  # now normalize across rows, by mean
  dmat <- dmat / apply(dmat, 1, mean, na.rm=T)
  # remove rows without quant
  dmat <- dmat[!apply(apply(dmat, 1, is.na), 2, sum) > 0,]
  dmat
}
normalize_ri_data_table <- function(ev, dcols, remove.empty.rows=T) {
  ev[,dcols][ev[,dcols]==0] <- NA
  
  ## normalize data
  
  # first normalize by column, by median
  ev[,dcols] <- t(t(ev[,dcols]) / apply(ev[,dcols], 2, median, na.rm=T))
  # now normalize across rows, by mean
  ev[,dcols] <- ev[,dcols] / apply(ev[,dcols], 1, mean, na.rm=T)
  # remove rows without quant
  if(remove.empty.rows) {
    ev <- ev[!apply(apply(ev[,dcols], 1, is.na), 2, sum) > 0,]
  }
  ev
}

extract_uniprot_id <- function(leading_protein) {
  sapply(strsplit(leading_protein, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })
}

# FROM APLPACK LIBRARY ------
faces<-function(xy,which.row,fill=FALSE,face.type=1,
                nrow.plot,ncol.plot,scale=TRUE,byrow=FALSE,main,
                labels,print.info = TRUE,na.rm = FALSE,
                ncolors=20,
                col.nose=rainbow(ncolors),                   # nose
                col.eyes=rainbow(ncolors,start=0.6,end=0.85),# eyes
                col.hair=terrain.colors(ncolors),            # hair
                col.face=heat.colors(ncolors),               # face
                col.lips=rainbow(ncolors,start=0.0,end=0.2), # lips
                col.ears=rainbow(ncolors,start=0.0,end=0.2), # ears
                
                plot.faces=TRUE, cex = 2){  # 180308 pwolf
  if((demo<-missing(xy))){
    xy<-rbind(
      c(1,3,5),c(3,5,7),
      c(1,5,3),c(3,7,5),
      c(3,1,5),c(5,3,7),
      c(3,5,1),c(5,7,3),
      c(5,1,3),c(7,3,5),
      c(5,3,1),c(7,5,3),
      c(1,1,1),c(4,4,4),c(7,7,7)
    )
    labels<-apply(xy,1,function(x) paste(x,collapse="-"))
  }  
  spline<-function(a,y,m=200,plot=FALSE){
    n<-length(a)
    h<-diff(a)
    dy<-diff(y)
    sigma<-dy/h
    lambda<-h[-1]/(hh<-h[-1]+h[-length(h)])
    mu<-1-lambda
    d<-6*diff(sigma)/hh
    tri.mat<-2*diag(n-2)
    tri.mat[2+  (0:(n-4))*(n-1)] <-mu[-1]
    tri.mat[    (1:(n-3))*(n-1)] <-lambda[-(n-2)]
    M<-c(0,solve(tri.mat)%*%d,0)
    x<-seq(from=a[1],to=a[n],length=m)
    anz.kl <- hist(x,breaks=a,plot=FALSE)$counts
    adj<-function(i) i-1
    i<-rep(1:(n-1),anz.kl)+1
    S.x<-  M[i-1]*(a[i]-x          )^3 / (6*h[adj(i)])  +
      M[i]  *(x        -a[i-1])^3 / (6*h[adj(i)])  +
      (y[i-1] - M[i-1]*h[adj(i)]^2 /6) * (a[i]-x)/ h[adj(i)] +
      (y[i]   - M[i]  *h[adj(i)]^2 /6) * (x-a[i-1]) / h[adj(i)]
    if(plot){ plot(x,S.x,type="l"); points(a,y)    }
    return(cbind(x,S.x))
  }
  
  n.char<-15
  xy<-rbind(xy)
  if(byrow) xy<-t(xy)
  if(any(is.na(xy))){
    if(na.rm){ 
      xy<-xy[!apply(is.na(xy),1,any),,drop=FALSE]
      if(nrow(xy)<3) {print("not enough data points"); return()}
      print("Warning: NA elements have been removed!!")
    }else{
      xy.means<-colMeans(xy,na.rm=TRUE)
      for(j in 1:length(xy[1,])) xy[is.na(xy[,j]),j]<-xy.means[j]
      print("Warning: NA elements have been exchanged by mean values!!")
    }  
  }
  if(!missing(which.row)&& all(  !is.na(match(which.row,1:dim(xy)[2]))  ))
    xy<-xy[,which.row,drop=FALSE]
  mm<-dim(xy)[2];  n<-dim(xy)[1]
  xnames<-dimnames(xy)[[1]]
  if(is.null(xnames)) xnames<-as.character(1:n)
  if(!missing(labels)) xnames<-labels
  if(scale){
    xy<-apply(xy,2,function(x){
      x<-x-min(x); x<-if(max(x)>0) 2*x/max(x)-1 else x })
  } else xy[]<-pmin(pmax(-1,xy),1)
  xy<-rbind(xy);n.c<-dim(xy)[2]
  # expand input matrix xy by replication of cols
  xy<-xy[,(rows.orig<-h<-rep(1:mm,ceiling(n.char/mm))),drop=FALSE]
  if(fill) xy[,-(1:n.c)]<-0
  
  face.orig<-list(
    eye  =rbind(c(12,0),c(19,8),c(30,8),c(37,0),c(30,-8),c(19,-8),c(12,0))
    ,iris =rbind(c(20,0),c(24,4),c(29,0),c(24,-5),c(20,0))
    ,lipso=rbind(c(0,-47),c( 7,-49),lipsiend=c( 16,-53),c( 7,-60),c(0,-62))
    ,lipsi=rbind(c(7,-54),c(0,-54))                  # add lipsiend
    ,nose =rbind(c(0,-6),c(3,-16),c(6,-30),c(0,-31))
    ,shape =rbind(c(0,44),c(29,40),c(51,22),hairend=c(54,11),earsta=c(52,-4),
                  earend=c(46,-36),c(38,-61),c(25,-83),c(0,-89))
    ,ear  =rbind(c(60,-11),c(57,-30))                # add earsta,earend
    ,hair =rbind(hair1=c(72,12),hair2=c(64,50),c(36,74),c(0,79)) # add hairend
  )
  lipso.refl.ind<-4:1
  lipsi.refl.ind<-1
  nose.refl.ind<-3:1
  hair.refl.ind<-3:1
  shape.refl.ind<-8:1
  shape.xnotnull<-2:8
  nose.xnotnull<-2:3
  
  nr<-n^0.5; nc<-n^0.5
  if(!missing(nrow.plot)) nr<-nrow.plot
  if(!missing(ncol.plot)) nc<-ncol.plot
  if(plot.faces){
    opar<-par(mfrow=c(ceiling(c(nr,nc))),oma=rep(6,4), mar=rep(.7,4))
    on.exit(par(opar))
  }
  
  face.list<-list()
  for(ind in 1:n){
    factors<-xy[ind,]
    face<-face.orig
    
    m<-mean(face$lipso[,2])
    face$lipso[,2]<-m+(face$lipso[,2]-m)*(1+0.7*factors[4])
    face$lipsi[,2]<-m+(face$lipsi[,2]-m)*(1+0.7*factors[4])
    face$lipso[,1]<-face$lipso[,1]*(1+0.7*factors[5])
    face$lipsi[,1]<-face$lipsi[,1]*(1+0.7*factors[5])
    face$lipso["lipsiend",2]<-face$lipso["lipsiend",2]+20*factors[6]
    
    m<-mean(face$eye[,2])
    face$eye[,2] <-m+(face$eye[,2] -m)*(1+0.7*factors[7])
    face$iris[,2]<-m+(face$iris[,2]-m)*(1+0.7*factors[7])
    m<-mean(face$eye[,1])
    face$eye[,1] <-m+(face$eye[,1] -m)*(1+0.7*factors[8])
    face$iris[,1]<-m+(face$iris[,1]-m)*(1+0.7*factors[8])
    
    
    m<-min(face$hair[,2])
    face$hair[,2]<-m+(face$hair[,2]-m)*(1+0.2*factors[9])
    m<-0
    face$hair[,1]<-m+(face$hair[,1]-m)*(1+0.2*factors[10])
    m<-0
    face$hair[c("hair1","hair2"),2]<-face$hair[c("hair1","hair2"),2]+50*factors[11]
    
    m<-mean(face$nose[,2])
    face$nose[,2]<-m+(face$nose[,2]-m)*(1+0.7*factors[12])
    face$nose[nose.xnotnull,1]<-face$nose[nose.xnotnull,1]*(1+factors[13])
    
    m<-mean(face$shape[c("earsta","earend"),1])
    face$ear[,1]<-m+(face$ear[,1]-m)* (1+0.7*factors[14])
    m<-min(face$ear[,2])
    face$ear[,2]<-m+(face$ear[,2]-m)* (1+0.7*factors[15])
    
    face<-lapply(face,function(x){ x[,2]<-x[,2]*(1+0.2*factors[1]);x})
    face<-lapply(face,function(x){ x[,1]<-x[,1]*(1+0.2*factors[2]);x})
    face<-lapply(face,function(x){ x[,1]<-ifelse(x[,1]>0,
                                                 ifelse(x[,2] > -30, x[,1],
                                                        pmax(0,x[,1]+(x[,2]+50)*0.2*sin(1.5*(-factors[3])))),0);x})
    #face$shape[,2]<-face$shape[,2]*(1+0.2*factors[1])
    #face$shape[,1]<-face$shape[,1]*(1+0.2*factors[2])
    #face$shape[,1]<-face$shape[,1]<-ifelse(face$shape[,1]>0,
    #   ifelse(face$shape[,2] > -30, face$shape[,1],
    #      pmax(0,face$shape[,1]+(face$shape[,2]+50)*0.2*sin(1.5*(-factors[3])))),0)
    
    invert<-function(x) cbind(-x[,1],x[,2])
    face.obj<-list(
      eyer=face$eye
      ,eyel=invert(face$eye)
      ,irisr=face$iris
      ,irisl=invert(face$iris)
      ,lipso=rbind(face$lipso,invert(face$lipso[lipso.refl.ind,]))
      ,lipsi=rbind(face$lipso["lipsiend",],face$lipsi,
                   invert(face$lipsi[lipsi.refl.ind,,drop=FALSE]),
                   invert(face$lipso["lipsiend",,drop=FALSE]))
      ,earr=rbind(face$shape["earsta",],face$ear,face$shape["earend",])
      ,earl=invert(rbind(face$shape["earsta",],face$ear,face$shape["earend",]))
      ,nose=rbind(face$nose,invert(face$nose[nose.refl.ind,]))
      ,hair=rbind(face$shape["hairend",],face$hair,invert(face$hair[hair.refl.ind,]),
                  invert(face$shape["hairend",,drop=FALSE]))
      ,shape=rbind(face$shape,invert(face$shape[shape.refl.ind,]))
    )
    face.list<-c(face.list,list(face.obj))
    
    if(plot.faces){
      plot(1,type="n",xlim=c(-105,105)*1.1, axes=FALSE,
           ylab="",ylim=c(-105,105)*1.3)
      title(xnames[ind], cex.main = cex, xpd = NA) #180308
      f<-1+(ncolors-1)*(factors+1)/2 # translate factors into color numbers
      xtrans<-function(x){x};  ytrans<-function(y){y}
      for(obj.ind in seq(face.obj)[c(10:11,1:9)]) {
        x <-face.obj[[obj.ind]][,1]; y<-face.obj[[obj.ind]][,2]
        xx<-spline(1:length(x),x,40,FALSE)[,2]
        yy<-spline(1:length(y),y,40,FALSE)[,2]
        if(plot.faces){ 
          lines(xx,yy)
          if(face.type>0){
            if(obj.ind==10) 
              polygon(xtrans(xx),ytrans(yy),col=col.hair[ceiling(mean(f[9:11]))],xpd=NA) # hair
            if(obj.ind==11){ 
              polygon(xtrans(xx),ytrans(yy),col=col.face[ceiling(mean(f[1:2 ]))],xpd=NA) # face
              if(face.type==2){
                # beard
                for(zzz in seq(hhh<-max(face.obj[[8]][,1]),-hhh,length=30)){
                  hrx<-rnorm(8,zzz,2); hry<-0:7*-3*rnorm(1,3)+abs(hrx)^2/150
                  hry<-min(face.obj[[9]][,2])+hry
                  lines(xtrans(hrx),ytrans(hry),lwd=5,col="#eeeeee",xpd=NA)
                }
                ind<-which.max(xx); wx<-xx[ind]; ind<-which.max(yy); wy<-yy[ind]
                # edge of hat
                wxh<-wx<-seq(-wx,wx,length=20); wyh<-wy<-wy-(wx-mean(wx))^2/250+runif(20)*3
                lines(xtrans(wxh),ytrans(wyh)); wx<-c(wx,rev(wx)); wy<-c(wy-10,rev(wy)+20)
                wmxy1<-wmxy0<-c(min(wx),min(wy)+20)
                wmxy2<-wmxy3<-c(runif(1,wmxy0[1],-wmxy0[1]), wy[1]+100)
                wmxy1[2]<-0.5*(wmxy0[2]+wmxy3[2]); wmxy2[1]<-0.5*(wmxy2[1]+wmxy0[1])
                npxy<-20; pxy<-seq(0,1,length=npxy)
                gew<-outer(pxy,0:3,"^")*outer(1-pxy,3:0,"^")*
                  matrix(c(1,3,3,1),npxy,4,byrow=TRUE)
                wxl<-wmxy0[1]*gew[,1]+wmxy1[1]*gew[,2]+wmxy2[1]*gew[,3]+wmxy3[1]*gew[,4]
                wyl<-wmxy0[2]*gew[,1]+wmxy1[2]*gew[,2]+wmxy2[2]*gew[,3]+wmxy3[2]*gew[,4]
                lines(xtrans(wxl),ytrans(wyl),col="green")
                wmxy1[1]<- wmxy0[1]<- -wmxy0[1]
                wmxy1[2]<-0.5*(wmxy0[2]+wmxy3[2]); wmxy2[1]<-0.5*(wmxy2[1]+wmxy0[1])
                wxr<-wmxy0[1]*gew[,1]+wmxy1[1]*gew[,2]+wmxy2[1]*gew[,3]+wmxy3[1]*gew[,4]
                wyr<-wmxy0[2]*gew[,1]+wmxy1[2]*gew[,2]+wmxy2[2]*gew[,3]+wmxy3[2]*gew[,4]
                points(xtrans(wmxy3[1]),ytrans(wmxy3[2]),pch=19,cex=2,col="#ffffff",xpd=NA)
                points(xtrans(wmxy3[1]),ytrans(wmxy3[2]),pch=11,cex=2.53,col="red",xpd=NA)
                polygon(xtrans(c(wxl,rev(wxr))),ytrans(c(wyl,rev(wyr))),col="red",xpd=NA) # hat
                polygon(xtrans(wx),ytrans(wy),col="#ffffff",xpd=NA) # edge of hat
              }
              
            }
            xx<-xtrans(xx); yy<-ytrans(yy)
            if(obj.ind %in% 1:2) polygon(xx,yy,col="#eeeeee") # eyes without iris
            if(obj.ind %in% 3:4) polygon(xx,yy,col=col.eyes[ceiling(mean(f[7:8 ]))],xpd=NA) # eyes:iris
            if(obj.ind %in% 9)   polygon(xx,yy,col=col.nose[ceiling(mean(f[12:13]))],xpd=NA)# nose
            if(obj.ind %in% 5:6) polygon(xx,yy,col=col.lips[ceiling(mean(f[1:3]))],xpd=NA)  # lips
            if(obj.ind %in% 7:8) polygon(xx,yy,col=col.ears[ceiling(mean(f[14:15]))],xpd=NA)# ears
            
          }
        }
      }
    }
    
  }
  
  if(plot.faces&&!missing(main)){
    par(opar);par(mfrow=c(1,1))
    mtext(main, 3, 3, TRUE, 0.5)
    title(main)
  }
  
  info<-c(
    "var1"="height of face   ",
    "var2"="width of face    ",
    "var3"="structure of face",
    "var4"="height of mouth  ",
    "var5"="width of mouth   ",
    "var6"="smiling          ",
    "var7"="height of eyes   ",
    "var8"="width of eyes    ",
    "var9"="height of hair   ",
    "var10"="width of hair   ",
    "var11"="style of hair   ",
    "var12"="height of nose  ",
    "var13"="width of nose   ",
    "var14"="width of ear    ",
    "var15"="height of ear   ")
  var.names<-dimnames(xy)[[2]]
  if(0==length(var.names)) var.names<-paste("Var",rows.orig,sep="")
  info<-cbind("modified item"=info,"Var"=var.names[1:length(info)])
  rownames(info)<-rep("",15)
  if(print.info){
    cat("effect of variables:\n")
    print(info)
  }
  if(demo&&plot.faces) {
    plot(1:15,1:15,type="n",axes=FALSE,bty="n")
    text(rep(1,15),15:1,adj=0,apply(info,1,function(x) 
      paste(x,collapse="   -   ")),cex=0.7)
  }
  
  names(face.list)<-xnames
  out<-list(faces=face.list,info=info,xy=t(xy))
  class(out)<-"faces"
  invisible(out)
  
}

