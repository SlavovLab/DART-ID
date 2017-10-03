library(readr)
#ev <- read_tsv('dat/ev.adj.elite.txt')
ev <- read_tsv('dat/evidence.txt')

exps <- unique(ev$`Raw file`)

cor.mat <- matrix(nrow=length(exps), ncol=length(exps))

for(i in exps) {
  ind.i <- match(i, exps)
  ev.i <- subset(ev, ev$`Raw file`==i & ev$PEP < 0.05)
  for(j in exps) {
    cat(((ind.i-1) * length(exps)) + ind.j, '/', length(exps) * length(exps), '\r')
    flush.console()
    
    ind.j <- match(j, exps)
    ev.j <- subset(ev, ev$`Raw file`==j & ev$PEP < 0.05)
    if(ind.i == ind.j) {
      # if the same, the correlation is 1
      cor.mat[ind.i, ind.j] <- 1
      next
    } 
      
    # get peptides shared between them
    shared.peps <- intersect(ev.i$Sequence, ev.j$Sequence)
    
    if(length(shared.peps) < 10) {
      cor.mat[ind.i, ind.j] <- NA
      next
    }
    
    rts <- matrix(unlist(lapply(shared.peps, FUN=function(sequence) {
      c(median(ev.i$`Retention time`[ev.i$Sequence==sequence]),
        median(ev.j$`Retention time`[ev.j$Sequence==sequence]))
    })), nrow=length(shared.peps), ncol=2, byrow=TRUE)
    
    cor.mat[ind.i, ind.j] <- cor(rts[,1], rts[,2])
  }
}
cor.mat[is.na(cor.mat)] <- 0

cor.mat.f <- as.data.frame(as.vector(cor.mat))
cor.mat.f$Row <- rep(seq(1,length(exps)), each=length(exps))
cor.mat.f$Col <- rep(seq(1,length(exps)), length(exps))
names(cor.mat.f)[1] <- 'Cor'

ggplot(cor.mat.f, aes(x=Col, y=Row, fill=Cor)) +
  geom_tile() +
  scale_y_continuous(trans='reverse')


hc.cor.mat <- hclust(dist(cor.mat))
exps.f <- exps[cutree(hc.cor.mat, k=10)==1]

image(cor.mat[hc.cor.mat$order, hc.cor.mat$order], col=heat.colors(12))

##
exps <- as.character(unlist(read_csv('dat/exps.corr.txt')[,2]))
cor.mat <- as.matrix(read_csv('dat/cor.mat.csv'))
cor.mat <- cor.mat[,-1]
rownames(cor.mat) <- NULL
colnames(cor.mat) <- NULL


