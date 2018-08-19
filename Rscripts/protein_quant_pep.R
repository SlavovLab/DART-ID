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

## filter data -----

dcols <- colnames(ev)[grepl('Reporter intensity corrected', colnames(ev))]

ev.f <- ev %>%
  # get UniProt accession ID
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))]) %>%
  # exclude CON, REV proteins
  filter(!grepl("CON__|REV__", Protein)) %>%
  # only select SQC master sets
  filter(grepl("SQC", `Raw file`)) %>%
  # filter out non J/U sets (QC40E,F, QC42A)
  filter(!grepl('SQC9', `Raw file`)) %>%
  # select 1% protein fdr
  #filter(!is.na(prot_fdr)) %>%
  #filter(prot_fdr < 0.01) %>%
  dplyr::select(c('Sequence', 'Modified sequence', 'Protein', 'Raw file', 'Retention time', 
                  'PEP', 'pep_new', 'pep_updated', 'qval', 'qval_updated', dcols))


ev.f <- ev.f %>%
  # remove empty channels
  dplyr::select(-grep('Reporter intensity corrected', colnames(ev.f))[c(3, 4)])

# only take rows w/ quantitation
ev.f <- ev.f[apply(ev.f[,grep('Reporter', colnames(ev.f))] == 0, 1, sum) == 0,]

# normalize data (by col, then by row)
dcols <- grep('Reporter intensity corrected', colnames(ev.f))
ev.f <- normalize_ri_data_table(ev.f, dcols)

prots <- unique(ev.f %>% filter(PEP < 1e-2) %>% pull(Protein))
prots <- intersect(prots, unique(ev.f %>% filter(PEP > 1e-2 & PEP < 3e-2) %>% pull(Protein)))

## filter, normalize, collapse, cluster data -------

ev_a <- ev.f %>%
  # only select from protein list
  filter(Protein %in% prots) %>%
  # filter at 1% PEP
  filter(PEP < 1e-2)

ev_b <- ev.f %>%
  filter(Protein %in% prots) %>%
  filter(PEP > 1e-2 & PEP < 3e-2)

ev_c <- ev.f %>% filter(PEP < 1e-2)
n <- sample.int(nrow(ev_c), size=floor(nrow(ev_c) / 2))
ev_d <- ev_c[n,]
ev_c <- ev_c[-n,]
cprots <- intersect(ev_c$Protein, ev_d$Protein)
ev_c <- ev_c %>% filter(Protein %in% cprots)
ev_d <- ev_d %>% filter(Protein %in% cprots)

dmat_a <- ev_a %>%
  # collapse data by sequence, by mean
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_a)[dcols], funs(mean)) %>%
  # collapse data by protein, by mean
  group_by(Protein) %>%
  summarise_at(colnames(ev_a)[dcols], funs(mean)) %>%
  # throw out the protein names
  dplyr::select(colnames(ev_a)[dcols])

dmat_b <- ev_b %>%
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_b)[dcols], funs(mean)) %>%
  group_by(Protein) %>%
  summarise_at(colnames(ev_b)[dcols], funs(mean)) %>%
  dplyr::select(colnames(ev_b)[dcols])

dmat_c <- ev_c %>%
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_c)[dcols], funs(mean)) %>%
  group_by(Protein) %>%
  summarise_at(colnames(ev_c)[dcols], funs(mean)) %>%
  dplyr::select(colnames(ev_c)[dcols])
dmat_d <- ev_d %>%
  group_by(`Modified sequence`, Protein) %>%
  summarise_at(colnames(ev_d)[dcols], funs(mean)) %>%
  group_by(Protein) %>%
  summarise_at(colnames(ev_d)[dcols], funs(mean)) %>%
  dplyr::select(colnames(ev_d)[dcols])

# cast to matrix, and manually cluster (J, then U)
dmat_a <- data.matrix(dmat_a)
dmat_b <- data.matrix(dmat_b)
dmat_c <- data.matrix(dmat_c)
dmat_d <- data.matrix(dmat_d)
dmat_a <- dmat_a[,c(1, 3, 5, 7, 2, 4, 6, 8)]
dmat_b <- dmat_b[,c(1, 3, 5, 7, 2, 4, 6, 8)]
dmat_c <- dmat_c[,c(1, 3, 5, 7, 2, 4, 6, 8)]
dmat_d <- dmat_d[,c(1, 3, 5, 7, 2, 4, 6, 8)]

j_channels=c(3, 5, 7)
u_channels=c(4, 6, 8)

# combine j and u channels to get 1 j/u vector for each dmat
ju_a <- apply(dmat_a[,j_channels], 1, mean) / apply(dmat_a[,u_channels], 1, mean)
ju_b <- apply(dmat_b[,j_channels], 1, mean) / apply(dmat_b[,u_channels], 1, mean)
ju_c <- apply(dmat_c[,j_channels], 1, mean) / apply(dmat_c[,u_channels], 1, mean)
ju_d <- apply(dmat_d[,j_channels], 1, mean) / apply(dmat_d[,u_channels], 1, mean)


# a -----------------------------------------------------------------------

pdf(file='Figures/ju_ratio_slavov.pdf', width=7, height=4)

layout(cbind(1, 2))

par(mar=c(3, 3, 0.1, 1),
    oma=c(0,0,2,0), pty='s')

plot(ju_a, ju_b, pch=16, cex=1, col=rgb(0,0,0,0.15),
     xlim=c(0.4, 2.2), ylim=c(0.4, 2.2), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
#     main='J/U Ratios, SQC Master Sets')
abline(a=0, b=1, col='red', lwd=2)
text(0.45, 2.1, bquote(rho*' = '*.(formatC(cor(ju_a, ju_b), digits=3))), 
     adj=c(0, 0.5), cex=1.2)
#text(1.45, 0.5, paste0(length(prots), ' proteins'), adj=c(0, 0.5), cex=1.2)

axis(1, at=seq(0, 2.5, by=0.5), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 2.5, by=0.5), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext(expression('PEP'<1e-2), side=1, line=1.5, cex=1)
mtext(expression(paste(1e-2<'PEP') < 3e-2), side=2, line=2, cex=1)

plot(ju_c, ju_d, pch=16, cex=1, col=rgb(0,0,1,0.15),
     xlim=c(0.4, 2.2), ylim=c(0.4, 2.2), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
abline(a=0, b=1, col='red', lwd=2)
text(0.45, 2.1, bquote(rho*' = '*.(formatC(cor(ju_c, ju_d), digits=3))), 
     adj=c(0, 0.5), col=rgb(0,0,1), cex=1.2)
#text(1.45, 0.5, paste0(length(cprots), ' proteins'), adj=c(0, 0.5), col=rgb(0,0,1), cex=1.2)

axis(1, at=seq(0, 2.5, by=0.5), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 2.5, by=0.5), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('PEP < 0.01 (Set 1)', side=1, line=1.5, cex=1)
mtext('PEP < 0.01 (Set 2)', side=2, line=2, cex=1)

mtext('Relative Protein Levels - Jurkat / U-937 Cells', side=3, 
      outer=T, cex=1.2, font=2, line=0.5)

dev.off()

# b -----------------------------------------------------------------------

# select only proteins that have high fold change in single cells
dmat_c <- cbind(dmat_a[,j_channels], dmat_b[,j_channels],
                dmat_a[,u_channels], dmat_b[,u_channels])
colnames(dmat_c) <- NA

# run t-test for each protein
prot_sigs <- sapply(1:length(prots), function(i) {
  t.test(dmat_c[i,1:(length(j_channels)*2)], 
         dmat_c[i,((length(j_channels)*2)+1):((length(j_channels)*2) + (length(u_channels)*2))], 
         alternative='two.sided', var.equal=T)$p.value
})
# apply bonferroni correction (FWER)
prot_sigs <- prot_sigs * length(prots)
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
png(file='manuscript/Figs/ratio_scatter_v3.png', width=7, height=6, units='in', res=250)

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
           col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))])
    
    abline(a=0, b=1, col='red', lwd=0.75)
    
    cor_text <- formatC(cor(log2(ratio_mat_a)[,i], log2(ratio_mat_b)[,j], method='pearson'), digits=3)
    text(-2.75, 2.5, bquote(.(as.name('rho'))*.('=')*.(cor_text)), adj=c(0, 0.5))
    
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
            side=3, cex=0.75, line=0.15)
    }
    if(j==4) {
      mtext(parse(text=paste0('J[',floor((i-1) / 2)+1,']/U[',((i-1) %% 2)+1,']')), 
            side=4, cex=0.75, line=0.5, las=1)
    }
  }
}

mtext('Log2 Jurkat / U937 Ratio − Spectra. PEP < 0.01', side=1, outer=T, cex=1, las=1, line=1.75)
mtext('Log2 Jurkat / U937 Ratio − DART-ID. PEP < 0.01', side=2, outer=T, cex=1, las=0, line=1.5)
mtext(paste0('Protein RI intensity, distinct peptides. n=', nrow(dmat_a2)), side=3, outer=T, cex=1, las=1, line=2, font=2)
#mtext('Top', side=3, outer=T, cex=1)
#mtext('Right', side=4, outer=T, cex=1, las=0)

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
