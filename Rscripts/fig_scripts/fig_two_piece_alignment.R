library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_1000cell_20180607_6/ev_updated.txt")

## two-piece alignment, examples -------

exps <- c(2, 3, 62)
betas <- c(4.390156572, 0.705147925, 0.79694858, 37.84193815, #2
           6.991690633, 0.864736612, 0.734386652, 26.73425284, #3
           -0.521905334, 0.998603863, 1.070621396, 37.30914361) #62
betas <- matrix(betas, nrow=3, byrow=T)

ev.f <- ev %>%
  filter(exp_id %in% exps) %>%
  filter(!is.na(muij))

seqs <- ev.f %>% 
  group_by(Sequence) %>% 
  summarise(n=length(unique(`Raw file`))) %>% 
  filter(n==length(exps)) %>%
  pull(Sequence)

#seqs <- seqs[1:75]
set.seed(1)
seqs <- seqs[sample.int(length(seqs), 75)]

ev.f <- ev.f %>%
  filter(Sequence %in% seqs)

## -------

pdf(file='manuscript/Figs/two_piece_alignment_v2.pdf', width=3.5, height=3)
par(mgp=c(1.2, 0.5, 0),
    mar=c(2.5,2,1.5,0.2),
    pin=c(2, 2),
    cex.lab=0.8, cex.axis=0.8, cex.main=1,
    #xaxs='i', yaxs='i',
    pty='s', las=1)

plot(0, 0, col='white', xlim=c(20, 40), ylim=c(20, 40),
     #xlab="Canonical RT (min)", ylab="Observed RT (min)",
     main='3 SCoPE-MS Experiments',
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n')

#pchs <- c(15, 16, 17)
pchs <- c('+', 'x', '*')
cexs <- c(1, 1, 2)
for(i in 1:length(exps)) {
  exp <- exps[i]
  
  segments(x0=0, x1=betas[i, 4], 
           y0=betas[i, 1], y1=betas[i, 1] + (betas[i, 2] * betas[i, 4]),
           col='red', lwd=1.5)
  segments(x0=betas[i, 4], x1=100, 
           y0=betas[i, 1] + (betas[i, 2] * betas[i, 4]), 
           y1=betas[i, 1] + (betas[i, 2] * betas[i, 4]) + ((100 - betas[i, 4]) * betas[i, 3]),
           col='blue', lwd=1.5)
  
  points(ev.f$mu[ev.f$exp_id==exp], ev.f$`Retention time`[ev.f$exp_id==exp], 
         pch=pchs[i], cex=cexs[i], col=rgb(0,0,0,0.4))
  
  #abline(v=betas[i, 4], col='green', lty=2, lwd=1.5)
}

text(x=28, y=21, labels='50 Shared Peptides', 
     cex=0.8, adj=c(0,1))

axis(side=1, tck=-0.02, padj=-0.6)
axis(side=2, tck=-0.02, padj=0.2)

mtext('Canonical RT (min)', 1, line=1.3)
mtext('Observed RT (min)', 2, line=1.6, las=3)

legend('topleft', c('First Segment', 'Second Segment'), 
       pch=c(NA, NA), col=c('red', 'blue'), lty=c(1, 1),
       lwd=c(2, 2), pt.cex=1.5, 
       bty='n', box.lwd=0,
       cex=0.8, y.intersp=1.1, inset=c(0.02, -0.02))

dev.off()

## observed vs. predicted --------

# pdf(file='manuscript/Figs/SFig_7.pdf', width=3, height=3)
# par(mgp=c(1.2, 0.5, 0),
#     mar=c(2.5,2,0.5,0.2),
#     pin=c(2, 2),
#     cex.lab=0.8, cex.axis=0.8, cex.main=0.8,
#     #xaxs='i', yaxs='i',
#     pty='s')
# 
# plot(0, 0, col='white', xlim=c(10, 50), ylim=c(10, 50),
#      xlab="Inferred RT (min)", ylab="Observed RT (min)",
#      xaxt='n', yaxt='n')
# for(i in 1:length(exps)) {
#   exp <- exps[i]
#   points(ev.f$muij[ev.f$exp_id==exp], ev.f$`Retention time`[ev.f$exp_id==exp], pch='x', cex=0.75)
# }
# abline(a=0, b=1, col='blue', lwd=1.5)
# 
# axis(side=1, tck=-0.02, padj=-0.6)
# axis(side=2, tck=-0.02, padj=0.2)
# 
# dev.off()
