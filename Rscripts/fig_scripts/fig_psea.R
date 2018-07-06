library(tidyverse)
source('Rscripts/lib.R')

# ev <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/ev_updated.txt')

## de-novo proteins ----------
# MSigDB


ev_a1 <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_v2.txt-1563089131_GO.txt', skip=0)
ev_b1 <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_denovo_v2.txt1471657536_GO.txt', skip=0)

ev_a2 <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_v2.txt-4388391_MSigDB.txt', skip=0)
ev_b2 <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_denovo_v2.txt-667074535_MSigDB.txt', skip=0)


ev_a1 <- ev_a1[-1,]
ev_b1 <- ev_b1[-1,]
ev_a2 <- ev_a2[-1,]
ev_b2 <- ev_b2[-1,]

ev_a1 <- ev_a1 %>% arrange(`Annotation Description`) %>% filter(`FDR (q-value)` < 0.25)
ev_b1 <- ev_b1 %>% arrange(`Annotation Description`) %>% filter(`FDR (q-value)` < 0.25)
ev_a2 <- ev_a2 %>% arrange(`Annotation Description`) %>% filter(`FDR` < 0.25)
ev_b2 <- ev_b2 %>% arrange(`Annotation Description`) %>% filter(`FDR` < 0.25)

common_annots_1 <- intersect(ev_a1$`Annotation Description`, ev_b1$`Annotation Description`)
common_annots_2 <- intersect(ev_a2$`Annotation Description`, ev_b2$`Annotation Description`)

ev_a1 <- ev_a1 %>% filter(`Annotation Description` %in% common_annots_1)
ev_b1 <- ev_b1 %>% filter(`Annotation Description` %in% common_annots_1)
ev_a2 <- ev_a2 %>% filter(`Annotation Description` %in% common_annots_2)
ev_b2 <- ev_b2 %>% filter(`Annotation Description` %in% common_annots_2)

## --------

pdf(file='manuscript/Figs/fig_psea.pdf', width=3, height=3)

par(mar=c(2,2.5,1.5,1),
    pty='s', las=1,
    cex.axis=0.75, cex.lab=0.75)

plot(0, 0, type='n',
     xlim=c(0, 6.75), ylim=c(0, 6.75),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')

# average # prots per group for both sets, then normalize to max

points(log2(ev_a1$PES), log2(ev_b1$PES), 
       pch=16, col=rgb(1, 0, 0, 0.5),
       cex=((ev_a1$`Number of proteins with annotation in dataset` + 
              ev_b1$`Number of proteins with annotation in dataset`) / 100) + 1)
points(log2(ev_a2$PES), log2(ev_b2$PES), 
       pch=16, col=rgb(0, 0, 1, 0.5),
       cex=((ev_a2$`Number of proteins with annotation in dataset` + 
              ev_b2$`Number of proteins with annotation in dataset`) / 100) + 1)
abline(a=0, b=1, col='black')

# legend('topleft', c('GO-Term', 'MSigDB', '', 'n=50', 'n=100', 'n=150'),
#        pch=16, col=c('red', 'blue', 'white', 'grey50','grey50','grey50'),
#        pt.cex=c(2,2,0,1.5, 2, 2.5), 
#        cex=0.8, ncol=2, x.intersp=1, y.intersp=1.4,
#        bty='n', inset=c(0.01, 0))
legend('topleft', c(paste0('GO-Term', ' | ', nrow(ev_a1), ' sets'),
                    paste0('MSigDB', ' | ', nrow(ev_a2), ' sets')), 
       pch=16, col=c('red', 'blue'), pt.cex=2, cex=0.8,
       bty='n', y.intersp=1.4, inset=c(0.02, 0))
legend('bottomright', c('n=50', 'n=100', 'n=150'),
       pch=16, col='grey50', pt.cex=c(1.5, 2, 2.5),
       cex=0.8, bty='n', y.intersp=1.4, inset=c(0, 0))

axis(1, at=seq(0, 7), tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 7), tck=-0.02, mgp=c(0, 0.4, 0))

mtext('Log2 Enrichment - Spectra', 1, line=1, cex=0.8)
mtext('Log2 Enrichment - DART-ID', 2, line=1, cex=0.8, las=3)
mtext('Protein Set Enrichment Analysis', 3, line=0.1, cex=0.8, font=2)

dev.off()

