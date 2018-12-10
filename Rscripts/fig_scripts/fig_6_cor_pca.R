library(tidyverse)
library(ggridges)
library(pracma)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## load cormat data -----

source('Rscripts/validation_cormats.R')

## ------

pdf(file='manuscript/Figs/cormats_type_v5.pdf', width=4.5, height=3)

# 1, 2 = cell type cor mats, 3 = colorbar
#layout(rbind(c(1, 1, 1, 3),
#             c(2, 2, 2, 3)))
layout(t(c(1, 2, 3)), widths=c(1, 1, 0.35))

# cell type cormats

colfunc <- colorRampPalette(c('blue', 'white', 'red'))
ncols <- 50

#cell_labels = parse(text=c('U937[3]', 'U937[2]', 'U937[1]',
#                           'Jurkat[3]', 'Jurkat[2]', 'Jurkat[1]'))
cell_labels=c('3', '2', '1', '3', '2', '1')

par(mar=c(0,0,0,0.5), 
    oma=c(1, 3, 3.5, 0),
    pty='s', cex.axis=1, cex.lab=1)

image(cor_type_a[nrow(cor_type_a):1,], zlim=c(-1, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')

axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)),
     labels=rev(cell_labels),
     tck=-0.02, las=0, mgp=c(0, 0.5, 0))
axis(2, at=seq(0, 1, length.out=nrow(cor_type_a)),
     labels=cell_labels,
     tck=-0.02, las=1, mgp=c(0, 0.5, 0))

mtext('Jurkat', side=1, at=0.2, line=1.75, cex=0.85, font=1)
mtext('U-937', side=1, at=0.8, line=1.75, cex=0.85, font=1)

mtext('Jurkat', side=2, at=0.8, line=1.5, cex=0.85)
mtext('U-937', side=2, at=0.2, line=1.5, cex=0.85)

#mtext(paste0('Spectra'), side=1, line=4, cex=1, font=2)

par(mar=c(0,0.5,0,0))
image(cor_type_b[nrow(cor_type_b):1,], zlim=c(-1, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')

axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)),
     labels=rev(cell_labels),
     tck=-0.02, las=0, mgp=c(0, 0.5, 0))
axis(2, at=seq(0, 1, length.out=nrow(cor_type_a)), labels=NA,
     tck=-0.02)

mtext('Jurkat', side=1, at=0.2, line=1.75, cex=0.85, font=1)
mtext('U-937', side=1, at=0.8, line=1.75, cex=0.85, font=1)

#mtext(paste0('DART-ID'), side=1, line=4, cex=1, font=2)

par(mar=c(2.75,0.75,2.75,2.5), pty='m')
image(matrix(seq(-1, 1, length.out=ncols), ncol=ncols), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(4, at=seq(0, 1, length.out=11), labels=seq(-1, 1, length.out=11), 
     tck=-0.1, las=1, mgp=c(0, 0.4, 0))

mtext('Pairwise correlations between cellular proteomes          ', 
      side=3, line=1, cex=1, font=2, las=1, outer=T)
mtext(paste0('Only proteins from Spectra\n', length(old_prots), ' proteins',
             ' | ', length(unique(ev_a$stan_peptide_id)), ' peptides'), 
      side=3, outer=T, line=-2.3, at=0.20, cex=0.85)
mtext(paste0('Only proteins from DART-ID\n', length(new_prots), ' proteins',
             ' | ', length(unique(ev_b$stan_peptide_id)), ' peptides'), 
      side=3, outer=T, line=-2.3, at=0.66, cex=0.85)

dev.off()


# distribution of correlations --------------------------------------------

cors_a <- as.vector(cor_type_a)
cors_a[cors_a == 1] <- NA
cors_b <- as.vector(cor_type_b)
cors_b[cors_b == 1] <- NA

df <- data.frame(cors=c(cors_a, cors_b),
                 method=rep(c('Spectra', 'DART-ID'), each=36))


# dot plot ----------------------------------------------------------------

p <- ggplot(df) +
  geom_dotplot(aes(x=method, y=cors, fill=method, color=method), 
               binaxis='y', stackdir='center', binwidth=0.05) +
  scale_x_discrete(limits=rev(levels(df$method))) +
  scale_y_continuous(limits=c(-0.75, 0.75), breaks=seq(-1, 1, by=0.25)) +
  scale_color_manual(values=c(cb[2], cb[1]), guide=F) +
  scale_fill_manual(values=c(cb[2], cb[1]), guide=F) +
  labs(x=NULL, y='Correlation') +
  theme_bert() + theme(
    plot.margin=margin(t=0.8,r=0.2,b=0.8,l=0.2,unit='cm'),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size=12),
    panel.grid.major.x = element_line(color=rgb(0,0,0,0.3), size=0.5),
    panel.grid.major.y = element_blank()
  )

ggsave(filename='manuscript/Figs/cortype_dotplot.pdf', plot=p, device='pdf', width=2.5, height=3)

# distribution plot --------------------------------------------------------------------

p <- ggplot(df) +
  #geom_density_ridges(aes(x=cors, y=method, group=method, fill=method), 
  #                    #stat='binline', bins=30,
  #                    rel_min_height=0.01, bandwidth=0.05) +
  geom_density_ridges_gradient(aes(x=cors, y=method, group=method, fill=..x..),
                               rel_min_height=0.01, bandwidth=0.05) +
  #scale_fill_manual(values=c(cb[2], cb[1]), guide=F) +
  scale_fill_gradient2(low='blue', mid='white', high='red', guide=F) +
  scale_x_continuous(limits=c(-1, 1)) +
  scale_y_discrete(expand=c(0.01, 0.01)) +
  labs(x='Correlation', y=NULL) +
  theme_ridges() + theme(
    plot.margin = margin(t=2, r=0.1, b=1.3, l=0.2, unit='cm'),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    axis.title.x = element_text(size=12, hjust=0.5, vjust=0.5)
  )
ggsave(filename='manuscript/Figs/cor_type_dists_v1.pdf', plot=p, device='pdf', width=2.5, height=3.5)

# ratios ------------------------------------------------------------------

pdf(file='manuscript/Figs/cormats_ratio.pdf', width=2, height=3)

# 1, 2 = j/u ratio cormats, 3 = colorbar
layout(rbind(c(1, 1, 1, 3),
             c(2, 2, 2, 3)))

par(oma=c(2, 0, 2, 0))

colfunc <- colorRampPalette(c('white', 'red'))
ratio_expressions <- c(
  'J[1]/U[1]', 'J[1]/U[2]', 'J[1]/U[3]',
  'J[2]/U[1]', 'J[2]/U[2]', 'J[2]/U[3]',
  'J[3]/U[1]', 'J[3]/U[2]', 'J[3]/U[3]')
ratio_expressions <- parse(text=ratio_expressions)

par(mar=c(0,3.5,0,0), pty='s',
    cex.axis=0.75, cex.lab=0.75)

image(cor_ratio_a[nrow(cor_ratio_a):1,], zlim=c(0.5, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     labels=NA, tck=-0.02)
axis(2, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=rev(ratio_expressions), mgp=c(0, 0.4, 0), las=1)

mtext('Log2 J/U Ratio Correlation', 3, line=0.5, cex=0.8, font=2, las=1)

par(mar=c(0,3.5,0,0))

image(cor_ratio_b[nrow(cor_ratio_b):1,], zlim=c(0.5, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=ratio_expressions, mgp=c(0, 0.3, 0), las=3)
axis(2, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=rev(ratio_expressions), mgp=c(0, 0.4, 0), las=1)

par(mar=c(2, 0.75, 2, 2), pty='m')
image(matrix(seq(-1, 1, length.out=ncols), ncol=ncols), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(4, at=seq(0, 1, length.out=6), labels=seq(0.5, 1, length.out=6), 
     tck=-0.1, las=1, mgp=c(0, 0.3, 0))

dev.off()



# Precursor Area, PIF dists -----------------------------------------------

ev_a <- ev.fi %>% filter(qval < 0.01)
ev_b <- ev.fi %>% filter(qval > 0.01 & qval_updated < 0.01)

# count missing data
dcols <- grep('corrected', colnames(ev.fi))[5:10]
dmat <- data.matrix(ev_a %>% select(colnames(ev.fi)[dcols]))
missing_a <- apply(dmat==0, 1, sum)
dmat <- data.matrix(ev_b %>% select(colnames(ev.fi)[dcols]))
missing_b <- apply(dmat==0, 1, sum)

# boxplots ----------------------------------------------------------------

pdf(file='manuscript/Figs/poor_quant.pdf', width=7, height=2)

layout(rbind(c(1, 2, 3, 4)))

par(cex.axis=1, 
    mar=c(3.5,2,1,2.5),
    oma=c(0,1,0,0))

boxplot(list(log10(ev_a$Intensity), log10(ev_b$Intensity)),
        range=1.5, col=c(cb[1], cb[2]), ylim=c(4.75, 8),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
#axis(1, at=c(1, 2), labels=c('Spectra', 'DART-ID'), tck=-0.02, mgp=c(0, 0.2, 0))
axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(c(1, 2), rep(4.3, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))

axis(2, at=seq(4, 10), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
mtext(expression('log'[10]*' Precursor Ion Area'), side=2, line=1.25, cex=0.85)

boxplot(list((ev_a$PIF), (ev_b$PIF)), 
        range=1.5, col=c(cb[1], cb[2]), ylim=c(0.55, 1),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0.3, 1, by=0.1), label=seq(0.3, 1, by=0.1)*100, tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2), rep(0.485, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Precursor Ion, %'), side=2, line=2, cex=0.85)

barplot(c(mean(ev_a$`Missed cleavages`), mean(ev_b$`Missed cleavages`)), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.2),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 0.25, by=0.05), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.02, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missed Cleavages, %'), side=2, line=2.25, cex=0.85)

barplot(c(mean(missing_a) / 6, mean(missing_b) / 6), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.15),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 1, by=0.05), labels=seq(0, 1, by=0.05)*100, 
     tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.015, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missing Data, %'), side=2, line=2, cex=0.85)

dev.off()


# PCA plots ---------

exps <- sort(unique(ev.f$`Raw file`))
# take 81 experiments 40-50, 80-90, 100-110
#exps <- exps[10:149]
# or exclude experiments with a lot of missing values
exps <- exps[-c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,22,27,28,33,38,39,54,55,56,57,60,61,62,
                64,67,68,69,76,77,78,80,81,94,95,102,103,108,110,111,112,113,114,115,116,
                117,120,121,122,123,124,125,126,127,128,129,130,132,133,137,139,147,149,150)]

# filter for the list of experiments
ev.fa <- ev.f %>% filter(`Raw file` %in% exps)

ev.fb <- ev.fa %>% filter(qval_updated < 0.01)
ev.fa <- ev.fa %>% filter(qval < 0.01)

prots.a <- ev.fa %>% 
  group_by(Protein) %>% 
  summarise(l=length(unique(`Raw file`))) %>% 
  filter(l >= length(exps)*0.95) %>%
  pull(Protein)

prots.b <- ev.fb %>% 
  group_by(Protein) %>% 
  summarise(l=length(unique(`Raw file`))) %>% 
  filter(l >= length(exps)*0.95) %>%
  pull(Protein)

prots.b <- prots.b[!prots.b %in% prots.a]

cat(length(prots.a), length(prots.b))

## -------

# select only proteins in protein list 
ev.fa <- ev.fa %>% filter(Protein %in% prots.a)
ev.fb <- ev.fb %>% filter(Protein %in% prots.b)

# number of peptides (only used in figure title)
peps.a <- unique(ev.fa$`Modified sequence`)
peps.b <- unique(ev.fb$`Modified sequence`)

dcols <- grep('Reporter intensity corrected', colnames(ev.fa))

dmat.fa <- ev.fa %>%
  # collapse data by sequence, by mean
  group_by(`Modified sequence`, Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fa)[dcols], funs(mean)) %>%
  # collapse data by protein, by mean
  group_by(Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fa)[dcols], funs(mean)) %>%
  ungroup()

dmat.fb <- ev.fb %>%
  group_by(`Modified sequence`, Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fb)[dcols], funs(mean)) %>%
  group_by(Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fb)[dcols], funs(mean)) %>%
  ungroup()

dmat.fa <- gather(dmat.fa, k, v, -c(Protein, `Raw file`))
dmat.fa <- spread(dmat.fa, 'Raw file', 'v')

dmat.fb <- gather(dmat.fb, k, v, -c(Protein, `Raw file`))
dmat.fb <- spread(dmat.fb, 'Raw file', 'v')

# impute
library(VIM)
dmat.fa[,-c(1:2)] <- VIM::kNN(dmat.fa[,-c(1:2)], imp_var=F)
dmat.fb[,-c(1:2)] <- VIM::kNN(dmat.fb[,-c(1:2)], imp_var=F)

## --------

# collapse experiments
cormat.a <- zeros(length(prots.a), length(exps)*length(dcols))
cormat.b <- zeros(length(prots.b), length(exps)*length(dcols))
for(i in 1:length(prots.a)) {
  cormat.a[i,] <- as.vector(data.matrix((dmat.fa %>% filter(Protein==prots.a[i]))[,-c(1,2)]))
}
for(i in 1:length(prots.b)) {
  cormat.b[i,] <- as.vector(data.matrix((dmat.fb %>% filter(Protein==prots.b[i]))[,-c(1,2)]))
}

# correlate cell types against each other
cormat.a <- cor(cormat.a)
cormat.b <- cor(cormat.b)

pars.a <- svd(cormat.a)
pars.b <- svd(cormat.b)



# plot PCA cells ----------------------------------------------------------

pdf(file='manuscript/Figs/pca.pdf', width=7, height=2.25)

par(oma=c(0, 2.5, 2.25, 2))

layout(rbind(c(1, 2)))

par(mar=c(2.4, 3.75, 0.5, 0.75), cex.axis=0.85)

plot(pars.a$u[,1], pars.a$u[,2], pch=16, cex=0.5,
     col=paste0(rep(c('#FF0000', '#0000FF'), nrow(pars.a$u)/2), '44'),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', xlim=c(-0.045, 0.045), ylim=c(-0.15, 0.17))

text(-0.019, 0.11, 'Jurkat', col='red', font=2, cex=1, adj=c(0, 0.5))
text(0.025, -0.1, 'U-937', col='blue', font=2, cex=1, adj=c(1, 0.5))

axis(1, at=seq(-0.1, 0.1, by=0.02), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(seq(-0.04, 0.04, by=0.02), -0.175, seq(-0.04, 0.04, by=0.02), 
     xpd=T, srt=30, adj=c(1, 1), cex=0.85)
axis(2, at=seq(-0.2, 0.2, by=0.05), las=1, tck=-0.02, mgp=c(0, 0.5, 0))

mtext(paste0('PC1 - ', formatC((pars.a$d[1] / sum(pars.a$d)) * 100, digits=3), '%'), side=1,
      cex=0.85, line=1.4)
mtext(paste0('PC2 - ', formatC((pars.a$d[2] / sum(pars.a$d)) * 100, digits=2), '%'), side=2,
      cex=0.85, line=2.5)
mtext('Only proteins from Spectra', side=3, cex=0.85, line=0.9)
mtext(paste0(length(prots.a), ' proteins | ', length(peps.a), ' peptides'), 
      side=3, cex=0.85, line=0.05)

#mtext(paste0('Separation of cellular proteomes - ',length(pars.a$d),' samples'), 
#      side=3, outer=T, cex=1, font=2, line=1.4)

mtext(paste0('Separation of cellular proteomes'), 
      side=3, outer=T, cex=1, font=2, line=1.4)


#par(mar=c(1.5, 1, 1, 1))

plot(pars.a$u[,1], pars.b$u[,2], pch=16, cex=0.5,
     col=paste0(rep(c('#FF0000', '#0000FF'), nrow(pars.b$u)/2), '44'),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', xlim=c(-0.045, 0.045), ylim=c(-0.15, 0.17))

text(-0.02, 0.11, 'Jurkat', col='red', font=2, cex=1, adj=c(0, 0.5))
text(0.025, -0.1, 'U-937', col='blue', font=2, cex=1, adj=c(1, 0.5))

axis(1, at=seq(-0.1, 0.1, by=0.02), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(seq(-0.04, 0.04, by=0.02), -0.175, seq(-0.04, 0.04, by=0.02), 
     xpd=T, srt=30, adj=c(1, 1), cex=0.85)
axis(2, at=seq(-0.2, 0.2, by=0.05), las=1, tck=-0.02, mgp=c(0, 0.5, 0))

mtext(paste0('PC1 - ', formatC((pars.b$d[1] / sum(pars.b$d)) * 100, digits=3), '%'), side=1,
      cex=0.85, line=1.4)
mtext(paste0('PC2 - ', formatC((pars.b$d[2] / sum(pars.b$d)) * 100, digits=2), '%'), side=2,
      cex=0.85, line=2.5)
mtext('Only proteins from DART-ID', side=3, cex=0.85, line=0.9)
mtext(paste0(length(prots.b), ' proteins | ', length(peps.b), ' peptides'), 
      side=3, cex=0.85, line=0.05)

dev.off()
