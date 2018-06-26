library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)

## Peptide Update ------

pdf(file='manuscript/Figs/confidence_update_v1.pdf', width=2.75, height=3)
par(mgp=c(1.2, 0.5, 0),
    mar=c(2.5,2,1.5,0.2),
    pin=c(2,2),
    cex.lab=0.8, cex.axis=0.8, cex.main=1,
    xaxs='i', yaxs='i',
    pty='s')

plot(0, 0,
     xlim=c(21, 25), ylim=c(-0.1, 1.85),
     xaxt='n', yaxt='n',
     xlab="Retention Time (min)", ylab="Density")
title(main="Confidence Update", line=0.5, adj=1)
#abline(v=mu1, col='blue', lty=2)
#segments(x0=mus[2]*exps[2], x1=mus[2]*exps[2], 
#         y0=0, y1=dnorm(mus[2]*exps[2], mean=mus[2]*exps[2], sd=0.275), col='blue', lty=2)
denx <- seq(21, 25, length.out=1000)
dist_sd <- 0.275
polygon(denx, dnorm(denx, mean=mus[9]*exps[3], sd=dist_sd), 
        col=rgb(0, 0, 1, 0.3))
lines(denx[400:1000], dnorm(denx[400:1000], mean=mus[9]*exps[3], sd=dist_sd), col='blue', lwd=2)

segments(x0=mus[9]*exps[3],x1=mus[9]*exps[3],
         y0=-1,y1=dnorm(mus[9]*exps[3], mean=mus[9]*exps[3], sd=dist_sd),
         col='blue', lwd=2, lty=2)
abline(h=0, col='black')

# null distribution
polygon(c(21,denx,25), c(0,(dlnorm(denx, meanlog=3.1, sdlog=0.07) * 1),0), 
        col=rgb(1, 0, 0, 0.3), border=NA)
lines(denx, dlnorm(denx, meanlog=3.1, sdlog=0.07) * 1, col='red', lwd=2)

xend=25.1

h2 <- 0.6
# segments(x0=x1,x1=x1,y0=0.07,y1=h2,
#          col='black', lwd=1.5)
# segments(x0=x1,x1=xend-0.4,y0=h2,y1=h2,
#          col='black', lwd=1.5)
# text(x=xend, y=h2, labels="PSM 1", adj=c(0.5, 0.5), cex=0.8)


h1 <- 0.45
# segments(x0=x2,x1=x2,y0=0.07,y1=h1-0.09,
#          col='black', lwd=1.5)
# #segments(x0=x2,x1=xend,y0=h1,y1=h1,
# #         col='black', lwd=1.5)
# text(x=xend,y=h1, labels="PSM 2", adj=c(0.5, 0.5), cex=0.8)

# points
x1 <- mus[9]*exps[3]-1.6
x2 <- mus[9]*exps[3]+0.3
points(c(x1, x2),rep(0,2), pch=c(16, 17), col='black', cex=1.5)

text(x=mus[9]*exps[3],y=1.65, labels="Peptide i\nExperiment A", font=1, cex=0.8)

axis(side=1, tck=-0.02, padj=-0.6)
axis(side=2, tck=-0.02, padj=0.2)

legend(x=21.1, y=0.5, xjust=0, yjust=0,
       c("Inferred RT", "Null RT", "Case 1", "Case 2"), pch=c(NA, NA, 16, 17),
       col=c("blue", "red", "black", "black"), lty=c(1, 1, NA, NA), lwd=c(2, 2, 5, 5),
       bty='n', cex=0.8, pt.cex=1.5, x.intersp=1, y.intersp=1.35,
       adj=c(0, 0.5), inset=c(0.1,0))

dev.off()

## bar plot -----

bdf <- data.frame(
  a=factor(c('Spectra', 'DART-ID'), levels=c('Spectra', 'DART-ID')),
  b=as.numeric(c(5e-2, 1e-3))
)

theme.bdf <- theme(
  axis.ticks.x=element_blank(),
  axis.title.y=element_text(size=10),
  axis.ticks.y=element_blank(),
  axis.text.y=element_blank(),
  plot.title=element_text(size=10, hjust=0.5, margin=margin(0,0,0.1,0,'cm')),
  legend.position=c(0.4,-0.2),
  legend.direction='horizontal',
  legend.key.size=unit(0.4,'cm'),
  legend.text=element_text(size=10),
  legend.margin=margin(0,0,0,0, unit='cm'),
  legend.spacing=unit(c(0,0,0,0), 'cm'),
  legend.box.spacing=unit(c(0,0,0,0), 'cm')
)

bar.plot.1 <- 
  ggplot(bdf) +
  geom_bar(aes(x=a, y=-log10(b), fill=a), stat='identity') +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.5), breaks=seq(0,3)) +
  scale_fill_manual(values=c(av[1],av[2]), guide=F) +
  labs(x=NULL, y='ID Confidence', fill=NULL, title="PSM 1") +
  #guides(fill=guide_legend(nrow=2)) +
  theme_bert() + theme.bdf + theme(
    plot.margin=margin(0.7,0.1,0.1,0.1,'cm')
  )

bdf$b <- c(5e-2, 0.3)
bar.plot.2 <- 
  ggplot(bdf) +
  geom_bar(aes(x=a, y=-log10(b), fill=a), stat='identity') +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.5), breaks=seq(0,3)) +
  scale_fill_manual(values=c(av[1],av[2]), guide=F) +
  labs(x=NULL, y='ID Confidence', fill=NULL, title="PSM 2") +
  guides(fill=guide_legend(nrow=2)) +
  theme_bert() + theme.bdf + theme(
    plot.margin=margin(0.1,0.1,1.35,0.1,'cm')
  )

bar.plot.1 <- ggplotGrob(bar.plot.1)
bar.plot.2 <- ggplotGrob(bar.plot.2)

gs <- list(bar.plot.1, bar.plot.2)
lay <- rbind(c(1),c(2))

pdf(file='manuscript/Figs/Fig_1D_v1.pdf', width=1.4, height=3.25)
grid.arrange(grobs=gs, layout_matrix=lay, heights=c(1, 1.2))
dev.off()

