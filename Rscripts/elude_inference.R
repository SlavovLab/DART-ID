library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

## load MaxQuant output ---------------------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ev_updated.txt')

# only take confident observations (from MaxQuant PEP, not DART-ID PEP)
# also remove SQC9*, which were run w/ a trapping column. The extra 10 or so minutes screws up
# some of the other prediction/alignment algos, so to be fair we are excluding those runs.
ev.f <- ev %>%
  filter(!is.na(pep_new)) %>%
  filter(!grepl('SQC9', `Raw file`))#%>%
  #filter(PEP < 0.01)
  
  
## load ssrcalc output ------------------------------------------------------------------------------
  hi_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA323_313_list1.txt')
  hi_dat <- rbind(hi_dat, read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA7794_648_list2.txt'))
  
  # append predicted HI to DART dataframe
  ev.f$HI <- hi_dat$`HI (pred)`[match(ev.f$Sequence, hi_dat$Sequence)]
  
  # global error
  ssrcalc_coefs <- lm(ev.f$`Retention time` ~ ev.f$HI)$coefficients
  ssrcalc_error <- ev.f$`Retention time` - ((ev.f$HI * ssrcalc_coefs[2]) + ssrcalc_coefs[1])
  
  ev.f$ssrcalc_error <- ssrcalc_error
  ev.f$ssrcalc_RT_corrected <- ((ev.f$HI * ssrcalc_coefs[2]) + ssrcalc_coefs[1])

## load ELUDE output ------------------------------------------------------------------------------

elude_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_out.txt', skip=2)
elude_dat <- elude_dat %>% 
  arrange(Peptide) %>% group_by(Peptide) %>%
  summarise(RT=unique(Predicted_RT))

# match up peptide sequences between MaxQuant output and ELUDE output
ev.f$elude_RT <- elude_dat$RT[match(ev.f$Sequence, elude_dat$Peptide)]

# remove NA observations, and run linear regression between observed RTs and predicted RTs,
# to correct for any global systematic biases in the LC configuration (such as added dwell time, etc.)
ev.a <- ev.f %>% filter(!is.na(elude_RT))
elude_coefs <- lm(ev.a$`Retention time` ~ ev.a$elude_RT)$coefficients
elude_error <- ev.f$`Retention time` - ((ev.f$elude_RT * elude_coefs[2]) + elude_coefs[1])

# apply linear regression coefs and put back into ev.f data frame
ev.f$elude_error <- elude_error
ev.f$elude_RT_corrected <- (ev.f$elude_RT * elude_coefs[2]) + elude_coefs[1]


# MaxQuant Match-Between-Runs (MBR) ---------------------------------------------------------------

exps <- sort(unique(ev.f$`Raw file`))
ev.mq <- data.frame()
for(exp in exps) {
  ev.a <- ev.f %>% 
    #filter(PEP < 0.01) %>% 
    filter(`Raw file`==exp) %>%
    filter(!is.na(`Calibrated retention time`))
  if(nrow(ev.a) == 0) return(c(0))
  coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
  ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.mq <- rbind(ev.mq, ev.a)
}

ev.a <- ev.f #%>% filter(PEP < 0.01)
coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])

mqmbr_error <- ev.mq$`Retention time`-ev.mq$RT_calibrated

# re-sort MBR dataframe by PSM id
ev.mq <- ev.mq %>% arrange(id)

ev.f$MBR_error <- mqmbr_error
ev.f$MBR_RT_corrected <- ev.mq$RT_calibrated

## refine PEPs with inferred RTs ------------------------------------------------------------------

# prune ev.f, only take the columns we need
ev.ff <- ev.f %>%
  select('Modified sequence', 'Raw file', 'Proteins', 'Leading razor protein', 'MS/MS scan number',
         'Retention time', 'PEP', 'pep_updated', 'pep_new',
         starts_with('Reporter intensity corrected'), 
         ends_with('RT_corrected'), ends_with('_error'))
  # remove rows w/o ELUDE data
  #filter(!is.na(elude_RT))

# null distribution: normal distribution with center at mean(RT) and standard deviation of sd(RT)
# rt_minus: probability of observing peptide RT at random. 
# i.e., density of RT evaluated on the null distribution
rt_minus <- dnorm(ev.ff$`Retention time`, mean=mean(ev.f$`Retention time`), sd=sd(ev.f$`Retention time`))

# rt_plus: probability of observing the peptide RT if the sequence is assigned correctly
# in this case, probability of observing RT with respect to ELUDE's predicted RT
# each unique peptide get's its own rt_plus distribution, where the mean is the ELUDE corrected RT,
# and the standard deviation is the standard deviation of the absolute ELUDE alignment error
rt_plus_elude <- dnorm(ev.ff$`Retention time`, mean=ev.ff$elude_RT_corrected, sd=sd(abs(ev.ff$elude_error), na.rm=T))
rt_plus_ssrcalc <- dnorm(ev.ff$`Retention time`, mean=ev.ff$ssrcalc_RT_corrected, sd=sd(abs(ev.ff$ssrcalc_error), na.rm=T))
rt_plus_MBR <- dnorm(ev.ff$`Retention time`, mean=ev.ff$MBR_RT_corrected, sd=sd(abs(ev.ff$MBR_error), na.rm=T))

# MaxQuant PEP estimates are strange and artefactual sometimes. Make sure none exceed 1 or are at 0
pep <- ev.ff$PEP
pep[pep > 1] <- 1
pep[pep == 0] <- 1e-100

#                                         P(RT|delta=0)*P(delta=0)
# PEP.new = P(delta=0|RT) =   ---------------------------------------------------
#                             P(RT|delta=0)*P(delta=0) + P(RT|delta=1)*P(delta=1)

elude_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_elude * (1 - pep)))
ev.ff$elude_pep <- elude_pep

ssrcalc_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_ssrcalc * (1 - pep)))
ev.ff$ssrcalc_pep <- ssrcalc_pep

MBR_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_MBR * (1 - pep)))
ev.ff$MBR_pep <- MBR_pep


## Define 2D PEP Scatter function, similar to Fig3A -----------------------------------------------


pep_scatter_plot <- function(pep, new_pep, filename, name, hi.color='black', title) {
  pep_log <- log10(pep)
  new_pep_log <- log10(new_pep)
  
  conf_limit <- 1e-8
  nbins <- 80
  x.bin <- seq(log10(conf_limit), 0, length=nbins)
  y.bin <- seq(log10(conf_limit), 0, length=nbins)
  
  freq <- as.data.frame(table(findInterval(pep_log, x.bin),
                              findInterval(new_pep_log, y.bin)))
  # remove any rows/columns over the 'nbins' limit
  freq <- freq %>%
    filter(as.numeric(Var1) <= nbins) %>%
    filter(as.numeric(Var2) <= nbins)
  
  freq[,1] <- as.numeric(freq[,1])
  freq[,2] <- as.numeric(freq[,2])
  
  freq2D <- diag(nbins)*0
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  
  colfunc <- colorRampPalette(c('white', hi.color))
  
  pdf(file=paste0('manuscript/Figs/',filename), width=2.5, height=2.5)
  
  layout(t(c(1, 2)), widths=c(6, 1))
  
  par(mar=c(1.5,2.5,1.25,0.75),
      pty='s', las=1,
      cex.axis=0.75, cex.lab=0.85, cex.main=1)
  
  cols <- colfunc(20)
  
  # Normal
  image(x.bin, y.bin, freq2D[-1,-1], col=cols,
        xlab=NA, ylab=NA,
        xaxs='i', yaxs='i',
        xaxt='n', yaxt='n', useRaster=F)
  
  abline(a=0, b=1, col='black')
  abline(h=-2, col='black', lty=2, lwd=1)
  abline(v=-2, col='black', lty=2, lwd=1)
  #segments(x0=-2, x1=-2, y0=log10(conf_limit), y1=-2, col='black', lty=2)
  #segments(x0=-2, x1=0, y0=-2, y1=-2, col='black', lty=2)
  
  rect(xleft=-2, xright=0, ybottom=log10(conf_limit), ytop=-2,
       border=NA, col=rgb(1,0,0,0.05))
  rect(xleft=log10(conf_limit), xright=-2, ybottom=-2, ytop=0,
       border=NA, col=rgb(0,0,1,0.05))
  
  text(-7, -1, 'Downgraded', cex=0.75, adj=c(0, 0.5))
  text(-1, -4.5, 'Upgraded', cex=0.75, adj=c(0, 0), srt=270)
  
  rng <- seq(-10, 0, 2)
  axis(1, tck=-0.02,  
       at=rng, labels=fancy_scientific(10^rng),
       mgp=c(0, 0.1, 0))
  axis(2, tck=-0.01, 
       at=rng, labels=fancy_scientific(10^rng),
       mgp=c(0, 0.3, 0), las=1)
  
  mtext('Spectra (MaxQuant)', 1, line=1.15, cex=0.85)
  mtext(name, 2, line=1.5, cex=0.85, las=3)
  mtext(paste0('Error Probability (PEP) - ', title), 3, line=0.6, cex=0.85, font=2)
  
  par(mar=c(2.5, 0.2, 3, 1), pty='m')
  image(matrix(seq(-1, 1, length.out=nbins), ncol=nbins), col=colfunc(nbins),
        xlab=NA, ylab=NA, xaxt='n', yaxt='n')
  mtext('Density', side=3, line=0.1, cex=0.75)
  
  dev.off()
}

# Generate plots

pep_scatter_plot(pep, ev.ff$pep_updated, 'pep_scatter_dart_v1.pdf', 'DART-ID', 'black', 'DART-ID')
pep_scatter_plot(pep, ev.ff$ssrcalc_pep, 'pep_scatter_ssrcalc_v1.pdf', 'SSRCalc', 'blue', 'SSRCalc')
pep_scatter_plot(pep, ev.ff$elude_pep, 'pep_scatter_elude_v1.pdf', 'ELUDE', 'red', 'ELUDE')
pep_scatter_plot(pep, ev.ff$MBR_pep, 'pep_scatter_MBR_v1.pdf', 'MaxQuant Calibrated RT', 'forestgreen', 'MaxQuant')

## Compare performance across these 4 methods -----------------------------------------------------

# first calculate q-values (FDR) from PEPs
pep_to_qval <- function(p) {
  (cumsum(p[order(p)]) / seq(1, length(p)))[order(order(p))]
}

ev.ff$qval_updated <- pep_to_qval(ev.ff$pep_updated)
ev.ff$ssrcalc_qval <- pep_to_qval(ev.ff$ssrcalc_pep)
ev.ff$elude_qval <- pep_to_qval(ev.ff$elude_pep)
ev.ff$MBR_qval <- pep_to_qval(ev.ff$MBR_pep)
ev.ff$qval <- pep_to_qval(pep)

# fold-change of IDs as a function of the confidence threshold
x <- logseq(5e-4, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'DART-ID', 'SSRCalc', 'ELUDE', 'MaxQuant MBR')
counter <- 1
for(i in x) {
  cat('\r', counter, '/', length(x), '       ')
  flush.console()
  counter <- counter + 1
  
  ratios <- c(
    1,
    sum(ev.ff$qval_updated < i) /    sum(ev.ff$qval < i),
    sum(ev.ff$ssrcalc_qval < i, na.rm=T) /    sum(ev.ff$qval < i & !is.na(ev.ff$ssrcalc_qval)),
    sum(ev.ff$elude_qval < i, na.rm=T)   /    sum(ev.ff$qval < i & !is.na(ev.ff$elude_qval)),
    sum(ev.ff$MBR_qval < i)     /    sum(ev.ff$qval < i  & !is.na(ev.ff$MBR_qval))
  )
  ident <- c(
    sum(ev.ff$qval < i) /              nrow(ev.ff),
    sum(ev.ff$qval_updated < i) /      nrow(ev.ff),
    sum(ev.ff$ssrcalc_qval < i, na.rm=T) /      sum(!is.na(ev.ff$ssrcalc_qval)),
    sum(ev.ff$elude_qval < i, na.rm=T) /        sum(!is.na(ev.ff$elude_qval)),
    sum(ev.ff$MBR_qval < i) /          sum(!is.na(ev.ff$MBR_qval))
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=method.names)


## FDR fold change --------------------------------------------------------------------------------

pdf(file='manuscript/Figs/inference_fdr_increase.pdf', width=2, height=4)

par(mar=c(2,2.75,2,0.15),
    #pty='s', 
    las=1,
    cex.axis=0.85, cex.lab=1, cex.main=1)

plot(0, 0, type='n',
     xlim=c(-3.1, -0.15), ylim=c(-15, 180),
     xlab=NA, ylab=NA, xaxs='i', yaxs='i', xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2, lwd=1)

cols <- c(rgb(0,0,0,0.5), 'black', 'blue', 'red', 'forestgreen')
#ltys <- c(1, 1, 2, 1, 2)
ltys <- c(1, 1, 1, 1, 1)
for(i in 1:length(method.names)) {
  df_a <- df %>% filter(Method == method.names[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lty=ltys[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.02,
     at=rng, 
     #labels=fancy_scientific(10^rng), 
     labels=c('0.1%', '1%', '10%', '100%'),
     mgp=c(0, -0.05, 0))
axis(2, tck=-0.01,
     #labels=NA,
     at=c(-25,seq(0, 175, 25)), 
     mgp=c(0, 0.2, 0))


legend('topright', c('Spectra', 'DART-ID', 'SSRCalc', 'ELUDE', 'MaxQuant\nMBR'),
       lwd=4, lty=1, col=c(rgb(0,0,0,0.5), 'black', 'blue', 'red', 'forestgreen'), seg.len=0.6,
       bty='n', cex=0.8, x.intersp=0.6, y.intersp=1.2, inset=c(-0.02, 0))

#mtext('FDR Threshold', 1, line=1, cex=1)

mtext('FDR Threshold', 1, line=1, cex=1)
mtext('% Increase', 2, line=1.5, cex=1, las=3)
mtext('Increase in PSMs', 3, line=0.2, cex=1, font=2)

dev.off()
