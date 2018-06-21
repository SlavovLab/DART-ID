## ---------

ev.f <- ev %>%
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  mutate_at('pep_updated', funs(ifelse(. > 1, 1, .))) %>%
  mutate(qval=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(ev)))[order(order(pep_updated))],
         dPEP=log10(PEP/pep_updated)) %>%
  filter(!grepl("CON__|REV__", Protein)) %>%
  filter(qval < 0.01) %>%
  filter(grepl("SQC", `Raw file`))

data.cols <- grep("Reporter intensity corrected", colnames(ev.f))
data.cols <- data.cols[-c(3,4)]

dmat <- data.matrix(ev.f[,data.cols])
# set zeroes to NA
dmat[dmat==0] <- NA

## normalize data

# first normalize by column, by median
dmat <- t(t(dmat) / apply(dmat, 2, median, na.rm=T))
# now normalize across rows, by mean
dmat <- dmat / apply(dmat, 1, mean, na.rm=T)
# remove rows without quant
dmat <- dmat[!apply(apply(dmat, 1, is.na), 2, sum) > 0,]

colnames(dmat) <- c('Carrier J', 'Carrier U', 'J', 'U', 'J', 'U', 'J', 'U')
heatmap(cor(dmat))

col_means <- as.numeric(apply(dmat, 2, mean))
sim <- dmat - matrix(rep(col_means, nrow(dmat)), nrow=nrow(dmat), byrow=T)
sim <- t(sim) %*% sim
heatmap(sim)

pars <- svd(sim)

par(xpd=T)
plot(pars$u[,1], pars$u[,2], pch=16,
     col=rep(c('red', 'blue'), 4), cex=c(3, 3, 1, 1, 1, 1, 1, 1),
     xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.75),
     xlab=paste0('PC1 - ', formatC(pars$d[1] / sum(pars$d), digits=3)),
     ylab=paste0('PC2 - ', formatC(pars$d[2] / sum(pars$d)), digits=3))
legend(x=0, y=0.75, c('J', 'U', 'Carrier'), pch=16, ncol=3,
       col=c('red', 'blue', 'black'), pt.cex=c(1, 1, 2), xjust=0.5, yjust=0)

pcamat <- pars$u
#colnames(pcamat) <- colnames(dmat)

faces(pcamat, labels=colnames(dmat), fill=F, face.type=1)
faces(cbind(pcamat,matrix(rep(1, 56), nrow=8)), labels=colnames(dmat), fill=T, face.type=1)
faces(sim)

## --------

facemat <- matrix(rep(1, 15*8), nrow=8)
facemat[,1:5] <- pcamat[,1:5]
facemat[,7:8] <- pcamat[,6:7]
facemat[,9] <- pcamat[,8]
facemat[,11] <- pcamat[,1]

faces(facemat, labels=colnames(dmat))
