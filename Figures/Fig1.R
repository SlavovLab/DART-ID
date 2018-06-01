## Figure 1 -------

library(tidyverse)
library(pracma)
library(gridExtra)
library(grid)
library(gtable)
library(RColorBrewer)
library(ggridges)
library(reshape2)
source('lib.R')

#ev <- read_tsv('dat/ev.adj.Fit3c.RTL.txt')
#ev <- read_tsv('dat/ev.adj.Fit3.txt')
#ev <- read_tsv('dat/ev.adj.Fit2.txt')

## Update Demo -------

f1.update.demo <- function() {

  x <- seq(0, 250, by=0.1)
  y1 <- dlnorm(x, meanlog=4, sdlog=0.5)
  y2 <- dnorm(x, mean=120, sd=2)
  
  #plot(x, y2, type='l', col='blue')
  #lines(x, y1, type='l', col='red')
  
  density.plot <- 
    ggplot(data.frame()) +
    geom_path(aes(x=x, y=y2, color='b'), size=0.75) +
    geom_path(aes(x=x, y=y1, color='a'), size=0.75) +
    scale_color_manual(values=c('red', 'blue', 'black')) +
    scale_x_continuous(limits=c(0, 170), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, 0.25), expand=c(0,0.003)) +
    labs(x=NULL, y='Density', color=NULL) +
    theme_bert() +
    theme(
      #plot.margin = unit(c(0.5,0.5,0.5,0.5), 'cm'),
      legend.position = c(0.1, 0.5),
      #legend.background = element_rect(fill='white', color=NULL, size=0.25),
      legend.key.height = unit(0.65, 'cm'),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  den.plot.1 <- 
    density.plot + 
    #geom_vline(xintercept=80, color='black', linetype='longdash') +
    geom_segment(aes(x=119, y=0, xend=119, yend=0.25, color='c'), size=0.5, linetype='longdash') +
    #annotate('text', x=80, y=-0.015, label='Observed RT', size=5) +
    scale_color_manual(values=c('red', 'blue', 'black'), guide=F)
  
  den.plot.2 <- 
    density.plot + 
    #geom_vline(xintercept=119, color='black', linetype='longdash') +
    geom_segment(aes(x=135, y=0, xend=135, yend=0.25, color='c'), size=0.5, linetype='longdash') +
    #annotate('text', x=119, y=-0.015, label='Observed RT', size=5) +
    scale_color_manual(values=c('red', 'blue', 'black'),
                       labels=c('Null RT Density',
                                'Predicted RT Density',
                                'Observed RT')) +
    guides(color=guide_legend(override.aes=list(linetype=c(1,1,2),
                                                size=c(0.75,0.75,0.5)))) +
    labs(y=' ') 
  den.plot.1 <- ggplotGrob(den.plot.1)
  den.plot.2 <- ggplotGrob(den.plot.2)
  
  
  bdf <- data.frame(
    a=factor(c('Spectra', 'Spectra+RT'), levels=c('Spectra', 'Spectra+RT')),
    b=as.numeric(c(5e-2, 1e-3))
  )
  
  bar.plot.1 <- 
  ggplot(bdf) +
    geom_bar(aes(x=a, y=-log10(b)), stat='identity',
             fill='grey90', color='black') +
    labs(x=NULL, y=parse(text='-log[10](PEP)')) +
    scale_y_continuous(expand=c(0,0), limits=c(0,3.5), breaks=seq(0,3)) +
    theme_bert() %+replace% theme(
      axis.text.x=element_text(size=11)
    )
  
  bdf$b <- c(5e-2, 0.3)
  bar.plot.2 <- 
    ggplot(bdf) +
    geom_bar(aes(x=a, y=-log10(b)), stat='identity',
             fill='grey90', color='black') +
    labs(x=NULL, y=' ') +
    scale_y_continuous(expand=c(0,0), limits=c(0,3.5), breaks=seq(0,3)) +
    theme_bert() %+replace% theme(
      axis.text.x=element_text(size=11),
      axis.line.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
  
  bar.plot.1 <- ggplotGrob(bar.plot.1)
  bar.plot.2 <- ggplotGrob(bar.plot.2)
  
  gs <- list(den.plot.1, den.plot.2, bar.plot.1, bar.plot.2,
             textGrob('       (i)', gp=gpar(fontsize=18)), textGrob('    (ii)', gp=gpar(fontsize=18)))
  lay <- rbind(c(1,2),c(3,4),c(5,6))
  
  #grid.arrange(den.plot.1, den.plot.2, bar.plot.1, bar.plot.2, 
  #             ncol=2, nrow=2, widths=c(1,1), heights=c(3,2))
  
  #pdf(file='manuscript/Figs/Fig_1B.pdf', width=5, height=4)
  
  grid.arrange(grobs=gs, layout_matrix=lay,
               widths=c(1,1), heights=c(3,2,0.5))
  
  #dev.off()
}


