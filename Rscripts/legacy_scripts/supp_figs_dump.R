## look at RT residuals of proline containing sequences ------

ev.f <- ev %>%
  filter(!is.na(muij)) %>%
  mutate(residual=`Retention time`-muij)

ev.f$proline <- "None"
ev.f$proline[grepl("P", ev.f$Sequence)] <- "Single"
ev.f$proline[grepl("PP", ev.f$Sequence)] <- "Double"
ev.f$proline[sample.int(nrow(ev.f), size=5000)] <- "Random"

sum(ev.f$proline=="Double")

ggplot(ev.f) +
  geom_density_ridges(aes(abs(residual),proline)) +
  scale_x_continuous(limits=c(0, 2.5)) +
  theme_ridges()

## SFig4 - resiudal analysis --------
# residual by RT bin

ev.f <- ev %>%
  filter(!is.na(muij)) %>%
  mutate(residual=`Retention time`-muij) %>%
  select(c("Sequence", "Raw file", "PEP", "pep_new", "pep_updated", 
           "Retention time", "muij", "residual")) %>%
  arrange(Sequence)

#ggplot(ev.f, aes(`Raw file`, residual)) +
#  geom_violin()


library(ggridges)
library(viridis)

ev.fa <- ev.f %>%
  mutate(rt_bin=cut(`Retention time`, breaks=c(0, 10, 20, 30, 40, 50, 60))) %>%
  filter(!is.na(rt_bin))

colfunc <- colorRampPalette(c('white', 'red'))

sfig4 <- 
  ggplot(ev.fa, aes(residual, rt_bin)) +
  geom_density_ridges(aes(fill=rt_bin), rel_min_height=0.01) +
  scale_x_continuous(limits=c(-1, 1), minor_breaks=seq(-1, 1, by=0.25)) +
  scale_fill_manual(values=colfunc(length(levels(ev.fa$rt_bin))), guide=F) +
  #scale_fill_manual(values=viridis(6), guide=F) +
  labs(x="Residual RT (min)", y="Retention time (min)") +
  theme_ridges() + theme(
    plot.margin=margin(0.2,0.1,0,0,'cm'),
    title=element_text(size=10),
    plot.title=element_text(size=10, face=1),
    panel.grid.minor=element_line(color=rgb(0,0,0,0.05)),
    axis.title.x=element_text(hjust=0.5),
    axis.title.y=element_text(hjust=0.5),
    axis.text=element_text(size=10)
  )
ggsave(plot=sfig4, file='manuscript/Figs/SFig_4.pdf', width=2, height=2.5)


## SFig5 -------

ev.fb <- ev.f %>% 
  filter(!is.na(residual)) %>%
  group_by(`Raw file`) %>% 
  summarise(res=median(abs(residual)),
            sigma=sd(residual-mean(residual)),
            res_05=quantile(abs(residual), 0.10),
            res_95=quantile(abs(residual), 0.90)) %>% 
  arrange(res) %>%
  mutate(x=seq(1,length(unique(ev.f$`Raw file`))),
         res_05=loess.smooth(seq(1, nrow(ev.fb)), res_05, evaluation=nrow(ev.fb))$y,
         res_95=loess.smooth(seq(1, nrow(ev.fb)), res_95, evaluation=nrow(ev.fb))$y)

sfig5 <-
  ggplot(ev.fb) +
  #geom_bar(aes(x, res), stat='identity')
  geom_ribbon(aes(x, ymin=res_05, ymax=res_95, fill="Confidence Interval, 10%-90%\n(Smoothed)")) +
  geom_path(aes(x, res, color="Mean |Residual RT|"), size=1) +
  #geom_smooth(aes(x, res_05), se=F) +
  #geom_smooth(aes(x, res_95), se=F) +
  #geom_path(aes(x, res_05)) +
  #geom_path(aes(x, res_95)) +
  scale_x_continuous(expand=c(0,0), breaks=seq(0, 200, by=30)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_manual(values=c(rgb(1,0,0,0.1))) +
  scale_color_manual(values=c('red'), labels=c("Median  | Residual RT |")) +
  guides(fill=guide_legend(order=2),
         color=guide_legend(order=1)) +
  labs(x="Experiment", y="| Residual RT | (min)",
       color=NULL, fill=NULL) +
  theme_bert() + theme(
    plot.margin=margin(0.2,0.1,0,0.1,'cm'),
    title=element_text(size=10),
    axis.text=element_text(size=10),
    plot.title=element_text(size=10, margin=margin(0,0,0,0,'cm'), hjust=0.5),
    legend.position=c(0.475, 0.9),
    legend.spacing.y=unit(0,'cm'),
    legend.margin=margin(0,0,0,0,'cm'),
    legend.key.height=unit(0.3,'cm'),
    legend.text=element_text(size=8)
  )
ggsave(plot=sfig5, file='manuscript/Figs/SFig_5.pdf', width=2.5, height=2.5)


## -------



## sfig3 ------------

pdf(file='manuscript/Figs/SFig_3.pdf', width=3, height=3)
par(mgp=c(1.2, 0.5, 0),
    mar=c(2.5,2,0.5,0.2),
    pin=c(2, 2),
    cex.lab=0.8, cex.axis=0.8, cex.main=0.8,
    #xaxs='i', yaxs='i',
    pty='s')

plot(df[df$x=="Spectra",]$y, df[df$x=="DART-ID",]$y,
     pch=16, col=paste0(av[2],"66"),
     xlab="# Spectral PSMs", ylab="# Updated PSMs", 
     #main="Peptides Quantified per Experiment",
     xlim=c(0, 1500), ylim=c(0, 2000),
     xaxt='n', yaxt='n')
points(df[df$x=="Spectra",]$y, df[df$x=="Percolator",]$y,
       pch=16, col=paste0(av[3],"66"))

abline(a=0, b=sum(df[df$x=="DART-ID",]$y) / sum(df[df$x=="Spectra",]$y), 
       col=paste0(av[2],"80"), lwd=2)
abline(a=0, b=sum(df[df$x=="Percolator",]$y) / sum(df[df$x=="Spectra",]$y), 
       col=paste0(av[3],"80"), lwd=2)

abline(a=0, b=1, col=av[1], lwd=3)

axis(side=1, at=seq(0, 1500, 250))
axis(side=2, at=seq(0, 2000, 500))
legend('bottomright', c("Spectra", "DART-ID", "Percolator"), 
       pch=c(NA, 16, 16), col=c(av[1], av[2], av[3]), lty=c(1, NA, NA),
       lwd=c(3, NA, NA), pt.cex=1.5, 
       bty='n', cex=0.8, y.intersp=1.1, inset=c(0.02, 0))
dev.off()

## --------

# df <- data.frame(
#   x=rep(apply(dmat, 2, sum), 2),
#   y=c(apply(dmat_new, 2, sum), apply(dmat_perc, 2, sum)),
#   Method=rep(c("Spectra+RT", "Percolator"), each=ncol(dmat))
# )
# ggplot(df, aes(x, y, color=Method)) +
#   geom_point() +
#   geom_abline(intercept=0, slope=1, color=av[1]) +
#   labs(color=NULL, x="Spectra", y="Updated") +
#   scale_color_manual(values=c(av[2], av[3])) +
#   guides(color=guide_legend(override.aes=list(size=4))) +
#   theme_bert() + theme(
#     aspect.ratio=1,
#     legend.position=c(0.7, 0.2),
#     legend.key.height=unit(0.5,'cm')
#   )
