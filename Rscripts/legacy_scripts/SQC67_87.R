library(tidyverse)
source('Rscripts/lib.R')
source('Figures/Fig1.R')
source('Figures/Fig2.R')

#ev <- read_tsv("~/git/RTLib/Alignments/SQC_1000cell_20180602_1/ev_updated.txt")
#ev <- read_tsv("~/git/RTLib/Alignments/SQC_1000cell_20180605_2/ev_updated.txt")
ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_1000cell_20180607_6/ev_updated.txt")

# remove excluded raw files to boost numbers
#ev <- ev %>% filter(!grepl("SQC72A|SQC72C|SQC72D|SQC74M|SQC67C16|SQC67C17|SQC67C18|SQC65|SQC95A9|SQC95A10|IFN6[H-K]-Trg|SQC67[AB]6", `Raw file`))

ev <- ev %>% filter(!grepl("IFN6[H-K]-Trg|SQC67[A-B]6|SQC67C1[3-9]|SQC72D", `Raw file`))

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

#fig2d <- f2.tmt.validation(ev)
fig2d <- f2.tmt.validation.3(ev)

fig2ef <- f2.protein.quant(ev)

gs <- list(ggplotGrob(fig2d), fig2ef[[1]], fig2ef[[2]])
lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3))
pdf(file='manuscript/Figs/Fig2DEF.pdf', width=7, height=2.75)
grid.arrange(grobs=gs, layout_matrix=lay)
dev.off()

## run with 1000 cells, and recent runs w/ trapping col

ev <- read_tsv("~/git/RTLib/Alignments/SQC_1000cell_20180602_1/ev_updated.txt")

ev.f <- ev %>% 
  filter(!is.na(pep_new)) %>%
  mutate(residual=muij-`Retention time`) %>%
  mutate(`With Trapping`=grepl("SQC95|SQC97", `Raw file`)) %>%
  mutate(`1000 cell`=grepl("FP18|IFN", `Raw file`))
  
ev.f %>% filter(abs(residual) < 2) %>% 
  #ggplot(aes(color=`With Trapping`)) +
  ggplot(aes(color=`1000 cell`)) +
  geom_freqpoly(aes(x=residual, y=..density..), size=1, position='identity', bins=100) +
  scale_x_continuous(limits=c(-2, 2)) + 
  labs(x="Residual RT (min)", y="Density",
       title=paste0("Residual RT\n", "n=", nrow(ev.f))) +
  guides(color=guide_legend(nrow=2)) +
  theme_bert() + theme(
    legend.position='bottom'
  )


## percolator analysis

pout <- read_tsv("SQCpout.tab")

# line up IDs
pout <- pout %>% arrange(PSMId)

#pout$PSMId %in% ev.f$id
ev <- cbind(ev, pout[match(ev$id, pout$PSMId),c("posterior_error_prob","q-value")])

ev <- ev %>%
  mutate(pep_perc=posterior_error_prob,
         fdr_perc=`q-value`)
ev$pep_perc_updated = ev$pep_perc
ev$pep_perc_updated[is.na(ev$pep_perc)] = ev$PEP[is.na(ev$pep_perc)]

ev %>% filter(!is.na(pep_perc)) %>% sample_n(1e4) %>%
ggplot() +
  geom_point(aes(x=PEP, y=pep_perc, color='Percolator PEP'), alpha=0.5, size=0.2) +
  geom_point(aes(x=PEP, y=pep_new, color='Spectra+RT PEP'), alpha=0.5, size=0.2)+
  geom_abline(aes(intercept=0, slope=1), color='black', size=1) +
  scale_x_log10(limits=c(1e-5, 1), expand=c(0,0)) +
  scale_y_log10(limits=c(1e-5, 1), expand=c(0,0)) +
  scale_color_manual(values=c('blue', 'red')) +
  guides(color=guide_legend(nrow=2, override.aes=list(size=5))) +
  labs(x="PEP", y="Updated PEP",
       title="PEP Update | RTLib vs. Percolator") +
  theme_bert() + theme(
    legend.position='bottom'
  )

ev %>% filter(!is.na(pep_perc)) %>% sample_n(1e4) %>%
  ggplot() +
  geom_point(aes(x=pep_perc, y=pep_new), color='blue', alpha=0.5, size=0.1) +
  geom_abline(aes(intercept=0, slope=1), color='black', size=1) +
  scale_x_log10(limits=c(1e-5, 1), expand=c(0,0)) +
  scale_y_log10(limits=c(1e-5, 1), expand=c(0,0)) +
  theme_bert()

## fig2abc -----

fig2a <- f2.pep.scatter(ev, bins=50)
ggsave(plot=fig2a, file='manuscript/Figs/Fig_2A_v1.pdf', width=2.32, height=2.25)

fig2bc <- f2.fold.change(ev)
fig2b <- fig2bc[[1]]
ggsave(plot=fig2b, file='manuscript/Figs/Fig_2B_v1.pdf', width=2.32, height=2.25)
fig2c <- fig2bc[[2]]
ggsave(plot=fig2c, file='manuscript/Figs/Fig_2C_v1.pdf', width=2.32, height=2.25)

# gs <- list(ggplotGrob(fig2a), ggplotGrob(fig2b), ggplotGrob(fig2c))
# lay <- rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
# pdf(file='manuscript/Figs/Fig2ABC.pdf', width=7, height=2.25)
# grid.arrange(grobs=gs, layout_matrix=lay)
# dev.off()

## fig2def -------

#fig2d <- ggplotGrob(f2.tmt.validation.3(ev))
#fig2d <- ggplotGrob(f2.tmt.validation.3(ev, cvs_all))
fig2d <- f2.tmt.validation.3(ev, cvs_all)
ggsave(plot=fig2d, file='manuscript/Figs/Fig_2D_v1.pdf', width=2.25, height=2.75)


fig2ef <- f2.protein.quant(ev)

ggsave(plot=fig2e, file='manuscript/Figs/Fig_2E_v1.pdf', width=3.5, height=2.75)
ggsave(plot=fig2f, file='manuscript/Figs/Fig_2F_v1.pdf', width=1.25, height=2.75)


#gs <- list(fig2d, ggplotGrob(fig2ef[[1]]), ggplotGrob(fig2ef[[2]]))
gs <- list(ggplotGrob(fig2d), ggplotGrob(fig2e), ggplotGrob(fig2f))
lay <- rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3))
pdf(file='manuscript/Figs/Fig2DEF.pdf', width=7, height=2.75)
grid.arrange(grobs=gs, layout_matrix=lay)
dev.off()



## FDR ------

ev$PEP[ev$PEP > 1] <- 1
ev$pep_updated[ev$pep_updated > 1] <- 1
#ev.f <- ev %>% filter(!is.na(pep_new))
fdr <- cumsum(ev$PEP[order(ev$PEP)]) / nrow(ev)
fdr_new <- cumsum(ev$pep_updated[order(ev$pep_updated)]) / nrow(ev)

x <- floor(seq(1, nrow(ev), length.out=1e4))
plot(1:length(x), fdr[x], type='l', col='blue', log='y',
     ylim=c(1e-10, 1))
lines(1:length(x), fdr_new[x], col='green')

## ------------

ev %>%
  mutate_at('pep_updated', funs(ifelse(. > 1, 1, .))) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev)))[order(order(PEP))],
         qval=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(ev)))[order(order(pep_updated))],
         dPEP=log10(PEP/pep_updated))


## ---------



