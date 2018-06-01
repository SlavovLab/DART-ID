library(tidyverse)
source('lib.R')
source('Figures/Fig1.R')
source('Figures/Fig2.R')

ev <- read_tsv("~/git/RTLib/Alignments/SQC67-87_20180531_2/ev_updated.txt")

ev <- ev %>% mutate(PEP.new=pep_new, PEP.updated=pep_updated)

#pdf(file='manuscript/Figs/Fig_1B.pdf', width=5, height=4)
fig1b <- f1.update.demo()
#dev.off()
ggsave('manuscript/Figs/Fig_1B.pdf', fig1b, device='pdf', width=5, height=4, units='in')

## -------
fig2a <- ggplotGrob(f2.pep.scatter(ev, bins=50))
#ggsave('manuscript/Figs/Fig_2A.pdf', fig2a, device='pdf', width=4, height=4, units='in')

fig2b <- ggplotGrob(f2.fold.change(ev))
#ggsave('manuscript/Figs/Fig_2B.pdf', fig2b, device='pdf', width=2, height=4, units='in')

fig2c <- ggplotGrob(f2.num.psms(ev))
#f2.num.psms(ev)
#ggsave('manuscript/Figs/Fig_2C.pdf', fig2c, device='pdf', width=2.5, height=4, units='in')

gs <- list(fig2a, fig2b, fig2c)
lay <- rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
pdf(file='manuscript/Figs/Fig2ABC.pdf', width=7, height=2.25)
grid.arrange(grobs=gs, layout_matrix=lay)
dev.off()

## ------

fig2d <- f2.tmt.validation(ev)

fig2ef <- f2.protein.quant(ev)

gs <- list(fig2d, fig2ef[1], fig2ef[2])
lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3))
pdf(file='manuscript/Figs/Fig2DEF.pdf', width=7, height=2.75)
grid.arrange(grobs=gs, layout_matrix=lay)
dev.off()
