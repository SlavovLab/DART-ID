library(tidyverse)

## de-novo proteins ----------
# MSigDB

ev_a <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_v2.txt-4388391_MSigDB.txt', skip=0)
#ev_a <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_v2.txt-1563089131_GO.txt', skip=0)
#ev_b <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_denovo_v2.txt-667074535_MSigDB.txt', skip=0)
#ev_b <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_denovo_v2.txt1471657536_GO.txt', skip=0)
ev_b <- read_tsv('~/git/RTLib/Alignments/SQC_20180621_2/PSEA_in_new_peptides.txt-624415452.txt')


ev_a <- ev_a[-1,]
ev_b <- ev_b[-1,]

ev_a <- ev_a %>% arrange(`Annotation Description`) %>% filter(`FDR` < 0.25)
ev_b <- ev_b %>% arrange(`Annotation Description`) %>% filter(`FDR` < 0.25)

common_annots <- intersect(ev_a$`Annotation Description`, ev_b$`Annotation Description`)

## --------

ev_a1 <- ev_a %>% filter(`Annotation Description` %in% common_annots)
ev_b1 <- ev_b %>% filter(`Annotation Description` %in% common_annots)

#avg_PES <- apply(cbind(ev_a1$PES, ev_b1$PES), 1, mean)
#palette <- colorRampPalette(c(''))

plot(ev_a1$PES, ev_b1$PES, pch=16, 
     main=paste('Protein Set Enrichment Analysis\nn=',length(common_annots)),
     xlab='PES - Proteins, FDR < 0.01',
     ylab='PES - New Proteins, Updated FDR < 0.01')
#text(ev_a1$PES+0.5, ev_b1$PES, labels=seq(1, length(common_annots)), cex=0.75, adj=c(0, 1))
abline(a=0, b=1, col='red')
