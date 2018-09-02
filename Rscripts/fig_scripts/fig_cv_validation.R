library(tidyverse)
library(ggridges)
source('Rscripts/lib.R')
source('Rscripts/validate.lib.3.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## add percolator data

source('Rscripts/add_percolator.R')

## run validation script -------

cvs_all <- validate.lib.3(ev, exclude.cols=c(1:4, 11))

#save(cvs_all, file='dat/cvs_all_20180712.rds')
#save(cvs_all, file='dat/cvs_all_20180813.rds')
save(cvs_all, file='dat/cvs_all_20180815.rds')


# load data ---------------------------------------------------------------

#load('dat/cvs_all_20180712.rds')
#load('dat/cvs_all_20180813.rds')
load('dat/cvs_all_20180815.rds')

## --------

methods <- as.character(unique(cvs_all$Method))
boxs <- list(DART=cvs_all$value[cvs_all$Method=='Spectra+RT'],
             Percolator=cvs_all$value[cvs_all$Method=='Percolator'],
             Spectra=cvs_all$value[cvs_all$Method=='Spectra'],
             Decoy=cvs_all$value[cvs_all$Method=='Null'])

# ridgeplot ---------------------------------------------------------------

#p <- 
ggplot(cvs_all) +
  geom_density_ridges(aes(x=value, y=rev(Method), group=rev(Method), fill=rev(Method)), 
                      rel_min_height=0.01, bandwidth=0.02) +
  scale_x_continuous(limits=c(0, 0.6)) +
  scale_y_discrete(limits=rev(levels(cvs_all$Method)), 
                   labels=rev(c('Decoy', 'Spectra', 'Percolator', 'DART-ID')),
                   expand=c(0.01, 0)) +
  scale_fill_manual(values=c(cb[4], cb[1], cb[3], cb[2]), guide=F) +
  labs(y=NULL, x='CV of Relative Quantification',
       title='Consistency of Protein Quantification') +
  theme_ridges() + theme(
    plot.margin=margin(0.25, 0.75, 0.25, 0.4, 'cm'),
    axis.title.x=element_text(hjust=0.5, size=12),
    #axis.title.y=element_text(hjust=0.5, size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    plot.title=element_text(size=12, hjust=0, vjust=1, lineheight=2, 
                            margin=margin(0,0,0.3,0,'cm'))
  )
ggsave('manuscript/Figs/cv_ridges_v2.pdf', p, 'pdf', width=4, height=2.25, units='in')

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

pdf(file='manuscript/Figs/cv_scatters_v4.pdf', width=4, height=2)

par(oma=c(0,1.25,0.5,0), 
    pty='s',
    cex.axis=0.85)

#layout(rbind(c(1), c(2)))
# layout(rbind(c(3, 7),
#              c(1, 4),
#              c(5, 7),
#              c(2, 6)),
#        heights=c(0.2, 1, 0.2, 1),
#        widths=c(1, 0.2))
layout(rbind(c(3, 7, 7, 5, 7),
             c(1, 4, 7, 2, 6)),
       heights=c(1, 5),
       widths=c(5, 1, 1, 5, 1))

par(mar=c(3.25, 2.25, 0.25, 0.25))

dens <- get_density(cvs_spectra$value, cvs_dart$value, k)

plot(cvs_spectra$value, cvs_dart$value, pch=16, cex=0.5,
     #col=rgb(0,0,0,0.2),
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
     xlim=c(0.05, 0.6), ylim=c(0.05, 0.6), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=0, b=1, col='red')

cor_text <- formatC(cor(cvs_spectra$value, cvs_dart$value), digits=3, format='f')
text(0.1, 0.55, bquote(rho*.(' = ')*.(cor_text)), adj=c(0, 0.5), cex=1)

axis(1, at=seq(0, 0.6, by=0.1), tck=-0.02, 
     #labels=NA,
     mgp=c(0, 0.1, 0))
axis(2, at=seq(0, 0.6, by=0.1), tck=-0.02,
     mgp=c(0, 0.4, 0), las=1)

mtext('Spectra CV', side=1, cex=0.85, line=1.4)
mtext('DART-ID CV', side=2, cex=0.85, line=1.7)

#par(mar=c(2.5, 2.25, 0.5, 0.25))
par(mar=c(3.25, 2.25, 0.25, 0.25), cex.axis=0.85)

dens <- get_density(cvs_spectra$value, cvs_decoy$value, k)

plot(cvs_spectra$value, cvs_decoy$value, pch=16, cex=0.5,
     #col=rgb(0,0,0,0.1),
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
     xlim=c(0.05, 0.6), ylim=c(0.05, 0.6), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=0, b=1, col='red')

cor_text <- formatC(cor(cvs_spectra$value, cvs_decoy$value), digits=3, format='f')
text(0.1, 0.55, bquote(rho*.(' = ')*.(cor_text)), adj=c(0, 0.5), cex=1)

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
par(mar=c(0,2.25,0.5,0.25), pty='m', xaxs='i', yaxs='i', cex.axis=0.85)
hist(cvs_spectra$value[cvs_spectra$value < 0.6], breaks=seq(0, 0.6, by=0.02), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
#axis(1, at=seq(0, 0.6, by=0.1), labels=NA, tck=-0.03)
axis(2, at=seq(0, 10, 3), labels=NA,
     tck=-0.03, mgp=c(0, 0.4, 0), las=1, cex=0.5)
mtext('Freq', side=2, cex=0.75, line=1.65)

# 4 = density of DART-ID CV
par(mar=c(3.25, 0, 0.25, 0.5), pty='m')
a <- hist(cvs_dart$value[cvs_dart$value < 0.6], breaks=seq(0, 0.6, by=0.02), plot=F)
barplot(a$density, space=0, horiz=T, axes=F, col='white')
axis(1, at=seq(0, 10, 3), labels=NA,
     tck=-0.03)
mtext('Freq', side=1, cex=0.75, line=1.25)

# 5 = density of spectra CV
par(mar=c(0,2.25,0.5,0.25), pty='m')
hist(cvs_spectra$value[cvs_spectra$value < 0.6], breaks=seq(0, 0.6, by=0.02), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
axis(2, at=seq(0, 10, 3), labels=NA,
     tck=-0.03, mgp=c(0, 0.4, 0), las=1, cex=0.5)
mtext('Freq', side=2, cex=0.75, line=1.65)

# 6 = density of decoy CV
par(mar=c(3.25, 0, 0.25, 0.5), pty='m')
a <- hist(cvs_decoy$value[cvs_decoy$value < 0.6], breaks=seq(0, 0.6, by=0.02), plot=F)
barplot(a$density, space=0, horiz=T, axes=F, col='white')
axis(1, at=seq(0, 10, 3), labels=NA,
     tck=-0.03)
mtext('Freq', side=1, cex=0.75, line=1.25)

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

