library(tidyverse)
source('Rscripts/lib.R')
source('Rscripts/validate.lib.3.R')

ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run validation script -------

cvs_all <- validate.lib.3(ev, exclude.cols=c(1:4, 11))
save(cvs_all, file='dat/cvs_all_20180712.rds')


# load data ---------------------------------------------------------------

load('dat/cvs_all_20180712.rds')

## --------

methods <- as.character(unique(cvs_all$Method))
boxs <- list(Spectra=cvs_all$value[cvs_all$Method=='Spectra'],
             Percolator=cvs_all$value[cvs_all$Method=='Percolator'],
             DART=cvs_all$value[cvs_all$Method=='Spectra+RT'],
             Decoy=cvs_all$value[cvs_all$Method=='Null'])

## --------

#pdf(file='manuscript/Figs/cv_validation.pdf', width=2, height=2.75)
pdf(file='manuscript/Figs/cv_validation_v3.pdf', width=1.5, height=5)

par(mar=c(4,3,0.25,0.5),
    oma=c(0, 0, 3, 0),
    pty='m', las=1,
    cex.axis=0.75, cex.lab=0.75)

boxplot(boxs, 
        col=c(av[1], av[3], av[2], av[5]),
        ylim=c(0.07, 0.4),
        xlab=NA, ylab=NA, xaxt='n', yaxt='n',
        outpch='x', outcex=1, outcol=rgb(0,0,0,0))

axis(1, at=seq(0, 6), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 1, 0.05), tck=-0.02, 
     mgp=c(0, 0.3, 0), las=1)
text(c(1, 2, 3, 4), par('usr')[3]-0.01, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID', 'Decoy'), cex=0.85)

mtext('CV of Relative Quantitation', 2, line=1.75, las=3, cex=1)
mtext('Consistency of\nProtein Quantitation', 3, line=0, cex=0.85, font=2, outer=T)

dev.off()


# scatter of cvs ----------------------------------------------------------

cvs_spectra <- cvs_all %>% filter(Method == 'Spectra')
cvs_dart <- cvs_all %>% filter(Method == 'Spectra+RT')
cvs_decoy <- cvs_all %>% filter(Method == 'Null')

k <- 80
contour_cols <- viridis::viridis(k, alpha=0.75)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

pdf(file='manuscript/Figs/cv_scatters_v2.pdf', width=2.5, height=5)

par(oma=c(0,1.25,0,0), 
    pty='s',
    cex.axis=0.75)

#layout(rbind(c(1), c(2)))
layout(rbind(c(3, 7),
             c(1, 4),
             c(5, 7),
             c(2, 6)),
       heights=c(0.2, 1, 0.2, 1),
       widths=c(1, 0.2))

par(mar=c(3.25, 2.25, 0.5, 0.1), cex.axis=0.85)

dens <- get_density(cvs_spectra$value, cvs_dart$value, k)

plot(cvs_spectra$value, cvs_dart$value, pch=16, cex=0.3,
     #col=rgb(0,0,0,0.2),
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
     xlim=c(0, 0.6), ylim=c(0, 0.6), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=0, b=1, col='red')

cor_text <- formatC(cor(cvs_spectra$value, cvs_dart$value), digits=3, format='f')
text(0.05, 0.75, bquote(rho*.(' = ')*.(cor_text)), adj=c(0, 0.5), cex=0.85)

axis(1, at=seq(0, 0.6, by=0.1), tck=-0.02, 
     #labels=NA,
     mgp=c(0, 0.1, 0))
axis(2, at=seq(0, 0.6, by=0.1), tck=-0.02,
     mgp=c(0, 0.4, 0), las=1)

mtext('Spectra CV', side=1, cex=0.85, line=1.4)
mtext('DART-ID CV', side=2, cex=0.85, line=1.7)

#par(mar=c(2.5, 2.25, 0.5, 0.25))
par(mar=c(3.25, 2.25, 0.5, 0.1), cex.axis=0.85)

dens <- get_density(cvs_spectra$value, cvs_decoy$value, k)

plot(cvs_spectra$value, cvs_decoy$value, pch=16, cex=0.3,
     #col=rgb(0,0,0,0.1),
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
     xlim=c(0, 0.6), ylim=c(0, 0.6), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=0, b=1, col='red')

cor_text <- formatC(cor(cvs_spectra$value, cvs_decoy$value), digits=3, format='f')
text(0.05, 0.75, bquote(rho*.(' = ')*.(cor_text)), adj=c(0, 0.5), cex=0.85)

axis(1, at=seq(0, 0.6, by=0.1), tck=-0.02,
     mgp=c(0, 0.1, 0))
axis(2, at=seq(0, 0.6, by=0.1), tck=-0.02,
     mgp=c(0, 0.4, 0), las=1)

mtext('Spectra CV', side=1, cex=0.85, line=1.4)
mtext('Decoy CV', side=2, cex=0.85, line=1.7)

#mtext('CV of Relative Quantitation', side=1, cex=1, line=0, outer=T)
#mtext('CV of Relative Quantitation', side=2, cex=1, line=1, outer=T)

# marginal densities
# 3 = density of spectra CV
par(mar=c(0,2.75,0.5,0.75), pty='m')
hist(cvs_spectra$value[cvs_spectra$value < 0.6], breaks=seq(0, 0.6, by=0.02), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
#axis(1, at=seq(0, 0.6, by=0.1), labels=NA, tck=-0.03)
axis(2, at=seq(0, 10, 3), labels=NA,
     tck=-0.03, mgp=c(0, 0.4, 0), las=1, cex=0.5)
mtext('Freq', side=2, cex=0.75, line=1.5)

# 4 = density of DART-ID CV
par(mar=c(3, 0, 0.25, 0.5), pty='m')
a <- hist(cvs_dart$value[cvs_dart$value < 0.6], breaks=seq(0, 0.6, by=0.02), plot=F)
barplot(a$density, space=0, horiz=T, axes=F, col='white')
axis(1, at=seq(0, 10, 3), labels=NA,
     tck=-0.03)
mtext('Freq', side=1, cex=0.75, line=1)

# 5 = density of spectra CV
par(mar=c(0,2.75,0.5,0.75), pty='m')
hist(cvs_spectra$value[cvs_spectra$value < 0.6], breaks=seq(0, 0.6, by=0.02), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
axis(2, at=seq(0, 10, 3), labels=NA,
     tck=-0.03, mgp=c(0, 0.4, 0), las=1, cex=0.5)
mtext('Freq', side=2, cex=0.75, line=1.5)

# 6 = density of decoy CV
par(mar=c(3, 0, 0.25, 0.5), pty='m')
a <- hist(cvs_decoy$value[cvs_decoy$value < 0.6], breaks=seq(0, 0.6, by=0.02), plot=F)
barplot(a$density, space=0, horiz=T, axes=F, col='white')
axis(1, at=seq(0, 10, 3), labels=NA,
     tck=-0.03)
mtext('Freq', side=1, cex=0.75, line=1)

dev.off()

## protein-by-protein --------

cvs_protein <- cvs_all %>% group_by(Var1, Method) %>%
  summarise(m=mean(value)) %>%
  spread(Method, m)

plot(cvs_protein$Spectra, cvs_protein$`Spectra+RT`, pch=16)
abline(a=0, b=1, col='red')

plot(cvs_protein$Spectra, cvs_protein$Percolator, pch=16)
abline(a=0, b=1, col='red')

plot(cvs_protein$Spectra, cvs_protein$Null, pch=16)
abline(a=0, b=1, col='red')


## --------

# t-tests
# t.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#        cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# t.test(cvs_all$value[cvs_all$Method=="DART-ID"], cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# t.test(cvs_all$value[cvs_all$Method=="Percolator"], cvs_all$value[cvs_all$Method=="Null"], var.equal=T)
# 
# var.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#          cvs_all$value[cvs_all$Method=="DART-ID"])
# t.test(cvs_all$value[cvs_all$Method=="Spectra"], 
#        cvs_all$value[cvs_all$Method=="DART-ID"], var.equal=T)

