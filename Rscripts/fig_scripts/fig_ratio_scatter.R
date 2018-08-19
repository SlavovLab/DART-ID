library(tidyverse)
library(MASS)
library(pracma)
library(RColorBrewer)
library(viridis)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180803_5exp_parametric_mixture_v2/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## load cormat data -----

source('Rscripts/validation_cormats.R')

## validation procedure ----------------------------
# need slightly different form of data. want to compare proteins now,
# don't want different proteins, but want different PSMs.

ev.f1 <- ev.f 
#ev.f1 <- normalize_ri_data_table(ev.f1, grep('Reporter intensity corrected', colnames(ev.f1)))

#ev_a <- ev.f %>% filter(qval < 0.001) %>%
ev_a <- ev.f1 %>% filter(PEP < 0.01)
#ev_b <- ev.f %>% filter(qval > 0.001 & qval_updated < 0.001) %>%
ev_b <- ev.f1 %>% filter(PEP > 0.01 & pep_updated < 0.01) %>%
  # select peptides that aren't in ev_a
  filter(!`Modified sequence` %in% unique(ev_a$`Modified sequence`))
  

dcols <- grep('Reporter intensity corrected', colnames(ev_a))
#ev_a <- normalize_ri_data_table(ev_a, dcols)
#ev_b <- normalize_ri_data_table(ev_b, dcols)

#common_peps <- intersect(ev_a$`Modified sequence`, ev_b$`Modified sequence`)
common_prots <- intersect(ev_a$Protein, ev_b$Protein)

collation_func <- mean

dmat_a <- ev_a %>%
  filter(Protein %in% common_prots) %>%
  #filter(`Modified sequence` %in% common_peps) %>%
  # collapse data by sequence, by mean
  group_by(`Modified sequence`, Protein) %>%
  #group_by(`Modified sequence`) %>%
  summarise_at(colnames(ev_a)[dcols], funs(collation_func)) %>%
  # collapse data by protein, by mean
  group_by(Protein) %>%
  summarise_at(colnames(ev_a)[dcols], funs(collation_func)) %>%
  # throw out the protein names
  dplyr::select(colnames(ev_a)[dcols])

dmat_b <- ev_b %>%
  filter(Protein %in% common_prots) %>%
  #filter(`Modified sequence` %in% common_peps) %>%
  group_by(`Modified sequence`, Protein) %>%
  #group_by(`Modified sequence`) %>%
  summarise_at(colnames(ev_b)[dcols], funs(collation_func)) %>%
  group_by(Protein) %>%
  summarise_at(colnames(ev_b)[dcols], funs(collation_func)) %>%
  dplyr::select(colnames(ev_b)[dcols])

# select only proteins that have high fold change in single cells
j_channels=c(3, 5, 7)
u_channels=c(4, 6, 8)

dmat_c <- cbind(dmat_a[,j_channels], dmat_b[,j_channels],
                dmat_a[,u_channels], dmat_b[,u_channels])
colnames(dmat_c) <- NA

# run t-test for each protein
prot_sigs <- sapply(1:length(common_prots), function(i) {
  t.test(dmat_c[i,1:(length(j_channels)*2)], 
         dmat_c[i,((length(j_channels)*2)+1):((length(j_channels)*2) + (length(u_channels)*2))], 
         alternative='two.sided', var.equal=T)$p.value
})
# apply bonferroni correction (FWER)
prot_sigs <- prot_sigs * length(common_prots)
fold_change_sig_thresh <- 0.05
cat(sum(prot_sigs < fold_change_sig_thresh), 'proteins with fold change significance p <', 
    fold_change_sig_thresh, 'after bonferroni correction')

dmat_a2 <- dmat_a[prot_sigs < fold_change_sig_thresh,]
dmat_b2 <- dmat_b[prot_sigs < fold_change_sig_thresh,]

# cast to matrix, and manually cluster (J, then U), dont select carrier channels
dmat_a2 <- data.matrix(dmat_a2)
dmat_b2 <- data.matrix(dmat_b2)
dmat_a2 <- dmat_a2[,c(3, 5, 7, 4, 6, 8)]
dmat_b2 <- dmat_b2[,c(3, 5, 7, 4, 6, 8)]

#j_channels <- c(2, 3, 4)
#u_channels <- c(6, 7, 8)
j_channels <- c(1, 2, 3)
u_channels <- c(4, 5, 6)

ratio_mat_a <- zeros(nrow(dmat_a2), length(j_channels) * length(u_channels))
ratio_mat_b <- zeros(nrow(dmat_b2), length(j_channels) * length(u_channels))
for(j in 1:length(j_channels)) {
  for(u in 1:length(u_channels)) {
    ratio_mat_a[,((j-1)*length(j_channels))+u] <- 
      dmat_a2[,j_channels[j]] / dmat_a2[,u_channels[u]]
    ratio_mat_b[,((j-1)*length(j_channels))+u] <- 
      dmat_b2[,j_channels[j]] / dmat_b2[,u_channels[u]]
  }
}

cor_ratio_a <- cor(log2(ratio_mat_a))
cor_ratio_b <- cor(log2(ratio_mat_b))
# cor_ratio_a <- cor(ratio_mat_a)
# cor_ratio_b <- cor(ratio_mat_b)

## log2 j/u ratio scatter -------

#pdf(file='manuscript/Figs/fig_ratio_scatter.pdf', width=5, height=5)
png(file='manuscript/Figs/ratio_scatter_v4.png', width=7, height=6, units='in', res=250)

par(mar=c(0.2,0.2,0.2,0.2),
    oma=c(3, 7, 4, 7), cex.axis=1,
    mfrow=c(4,4), pty='s')

k <- 60
contour_cols <- viridis(k)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

for(i in 1:4) {
  for(j in 1:4) {
    
    plot(0,0,type='n',
         xlim=c(-3, 3), ylim=c(-3, 3),
         xaxt='n', yaxt='n', xlab=NA, ylab=NA)
    
    abline(h=0, lwd=0.75, lty=2)
    abline(v=0, lwd=0.75, lty=2)
    
    dens <- get_density(log2(ratio_mat_a)[,i], log2(ratio_mat_b)[,j], k)
    
    points(log2(ratio_mat_a)[,i], log2(ratio_mat_b)[,j],
         pch=16, cex=0.75, 
         #col=rgb(0,0,0,0.2),
         #col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))]
         col=rgb(0,0,0,0.5)
         )
    
    abline(a=0, b=1, col='red', lwd=0.75)
    
    cor_text <- formatC(cor(log2(ratio_mat_a)[,i], log2(ratio_mat_b)[,j], method='pearson'), digits=3)
    text(-2.95, 2.5, bquote(.(as.name('rho'))*.(' = ')*.(cor_text)), adj=c(0, 0.5))
    
    # draw axes for marginal squares
    if(i == 4) {
      par(mgp=c(0, 0.2, 0))
      axis(1, at=c(-2, 0, 2), tck=-0.02)
    }
    if(j == 1) {
      par(mgp=c(0, 0.4, 0))
      axis(2, at=c(-2, 0, 2), tck=-0.02, las=1)
    }
    # label rows and columns on marginal squares
    if(i == 1) {
      mtext(parse(text=paste0('J[',floor((j-1) / 2)+1,']/U[',((j-1) %% 2)+1,']')), 
            side=3, cex=1, line=0.15)
    }
    if(j==4) {
      mtext(parse(text=paste0('J[',floor((i-1) / 2)+1,']/U[',((i-1) %% 2)+1,']')), 
            side=4, cex=1, line=0.5, las=1)
    }
  }
}

mtext('Spectral PEP < 0.01', side=1, outer=T, cex=1, las=1, line=1.75)
mtext('Spectral PEP > 0.01  &  DART-ID PEP < 0.01', side=2, outer=T, cex=1, las=0, line=1.5)
mtext(paste0('Log2 Jurkat / U-937 Reporter Ion Intensities'), side=3, 
      outer=T, cex=1, las=1, line=2, font=2)

dev.off()

## sampling to create decoy proteins/null distributions ------





## ---------

conf_thresh <- 0.01
ev_a <- ev.f %>%
  mutate(confident=as.numeric(PEP < conf_thresh),
         new_confident=as.numeric((PEP > conf_thresh) & (pep_new < conf_thresh))) %>%
  filter(confident > 0 | new_confident > 0) %>%
  dplyr::select(-grep('Reporter intensity corrected', colnames(ev.f))[c(3, 4)])

common_peps <- ev_a %>%
  group_by(`Modified sequence`) %>%
  summarise(n = any(confident == 1) & any(new_confident == 1),
            m = sum(confident),
            o = sum(new_confident)) %>%
  filter(n) %>%
  pull(`Modified sequence`)
common_prots <- ev_a %>%
  group_by(Protein) %>%
  summarise(n = (any(confident == 1) & any(new_confident == 1)),
            m = sum(confident),
            o = sum(new_confident)) %>%
  filter(n & m > 5 & o > 5) %>%
  pull(Protein)

ev_a <- ev_a %>% 
  filter(Protein %in% common_prots) %>%
  mutate(confident=confident + (new_confident*2)) %>%
  dplyr::select(-new_confident)

dcols <- grep('Reporter intensity corrected', colnames(ev_a))
ev_a <- normalize_ri_data_table(ev_a, dcols)

j_channels <- dcols[c(2, 3, 4)]
u_channels <- dcols[c(6, 7, 8)]

ratio_mat <- zeros(nrow(ev_a), length(j_channels) * length(u_channels))
for(j in 1:length(j_channels)) {
  for(u in 1:length(u_channels)) {
    ratio_mat[,((j-1)*length(j_channels))+u] <- 
      data.matrix(ev_a[,j_channels[j]]) / data.matrix(ev_a[,u_channels[u]])
  }
}
colnames(ratio_mat) <- paste0('ratio_', seq(1,length(j_channels)*length(u_channels)))
ev_a <- cbind(ev_a, ratio_mat)

ratio_cols <- grep('ratio', colnames(ev_a))
a <- ev_a %>%
  group_by(`Protein`, confident) %>%
  summarise_at(colnames(ev_a)[ratio_cols], funs(mean, sd))
b <- a %>% filter(confident == 2)
a <- a %>% filter(confident == 1)

plot(a$ratio_2_mean, b$ratio_2_mean)
abline(a=0, b=1)
