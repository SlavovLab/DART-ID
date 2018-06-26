library(tidyverse)
source('Rscripts/lib.R')
source('Rscripts/validate.lib.3.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run validation script -------

cvs_all <- validate.lib.3(ev, exclude.cols=c(1:4, 11))

## --------

methods <- as.character(unique(cvs_all$Method))
boxs <- list(Spectra=cvs_all$value[cvs_all$Method=='Spectra'],
             Percolator=cvs_all$value[cvs_all$Method=='Percolator'],
             DART=cvs_all$value[cvs_all$Method=='Spectra+RT'],
             Decoy=cvs_all$value[cvs_all$Method=='Null'])

## --------

pdf(file='manuscript/Figs/cv_validation.pdf', width=2, height=2.75)

par(mar=c(2.5,2.5,1.75,0.25),
    pty='m', las=1,
    cex.axis=0.75, cex.lab=0.75)

boxplot(boxs, 
        col=c(av[1], av[3], av[2], av[5]),
        ylim=c(0.1, 0.4),
        xlab=NA, ylab=NA, xaxt='n', yaxt='n',
        outpch='x', outcex=1, outcol=rgb(0,0,0,0))

axis(1, at=c(1, 2, 3, 4), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 1, 0.1), tck=-0.02, 
     mgp=c(0, 0.3, 0), las=1)
text(c(1, 2, 3, 4), par('usr')[3]-0.015, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID', 'Decoy'), cex=0.7)

mtext('CV of Relative Quantitation', 2, line=1.2, las=3, cex=0.75)
mtext('Consistency of\nProtein Quantitation', 3, line=0, cex=0.75, font=2)

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

