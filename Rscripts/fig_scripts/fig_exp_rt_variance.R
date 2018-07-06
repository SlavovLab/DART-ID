library(ggridges)
library(viridis)

ev.fa <- ev %>%
  filter(!is.na(muij)) %>%
  mutate(residual=`Retention time`-muij) %>%
  #mutate(rt_bin=cut(`Retention time`, breaks=c(0, 10, 20, 30, 40, 50, 60))) %>%
  mutate(rt_bin=cut(`Retention time`, breaks=seq(0, 60, by=10))) %>%
  filter(!is.na(rt_bin))

bin_density <- data.frame(table(ev.fa$rt_bin))
bin_density$Freq <- bin_density$Freq / sum(bin_density$Freq)

colfunc <- colorRampPalette(c('white', 'red'))

ggsave(ggplot(ev.fa, aes(residual, rt_bin)) +
  geom_density_ridges(aes(fill=bin_density$Freq[as.numeric(rt_bin)]), rel_min_height=0.01) +
  scale_x_continuous(limits=c(-1, 1), minor_breaks=seq(-1, 1, by=0.25)) +
  #scale_fill_manual(values=colfunc(length(levels(ev.fa$rt_bin))), guide=F) +
  scale_fill_continuous(low='white', high='red', limits=c(0, 0.4)) +
  labs(x="Residual RT (min)", y="Retention time (min)",
       fill='Fraction of\nall Data\n',
       title='Residual RT Increases With Time') +
  theme_ridges() + theme(
    plot.margin=margin(0.2,0.1,0,0,'cm'),
    title=element_text(size=10, hjust=0.5),
    plot.title=element_text(size=12, face='bold', 
                            margin=margin(0,0,0.1,-0.2,'cm'), hjust=0),
    panel.grid.minor=element_line(color=rgb(0,0,0,0.05)),
    #legend.margin =margin(0.5,0,0,0,'cm'),
    legend.title = element_text(margin=margin(0,0,5,0,'cm'), hjust=0.5),
    #legend.box.margin = margin(0.5,0,0,0,'cm'),
    #legend.box.spacing=unit(c(0.5,0,0,0),'cm'),
    #legend.margin = margin(0.5,0,0,0,'cm'),
    legend.text = element_text(size=10),
    axis.title.x=element_text(hjust=0.5),
    axis.title.y=element_text(hjust=0.5),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12)
  ), file='manuscript/Figs/fig_residual_rt_time.pdf', width=3.5, height=3)

## SFig5 -------

ev.fb <- ev %>% 
  filter(!is.na(muij)) %>%
  mutate(residual=`Retention time`-muij) %>%
  filter(!is.na(residual)) %>%
  group_by(`Raw file`) %>% 
  summarise(res=median(abs(residual)),
            sigma=sd(residual-mean(residual)),
            res_10=quantile(abs(residual), 0.10),
            res_90=quantile(abs(residual), 0.90),
            res_20=quantile(abs(residual), 0.20),
            res_80=quantile(abs(residual), 0.80)) %>% 
  arrange(res) %>%
  mutate(x=seq(1,length(unique(`Raw file`))),
         res_10=loess.smooth(seq(1, length(unique(`Raw file`))), res_10, 
                             evaluation=length(unique(`Raw file`)))$y,
         res_90=loess.smooth(seq(1, length(unique(`Raw file`))), res_90, 
                             evaluation=length(unique(`Raw file`)))$y,
         res_20=loess.smooth(seq(1, length(unique(`Raw file`))), res_20, 
                             evaluation=length(unique(`Raw file`)))$y,
         res_80=loess.smooth(seq(1, length(unique(`Raw file`))), res_80, 
                             evaluation=length(unique(`Raw file`)))$y)

## --------

pdf(file='manuscript/Figs/fig_residual_rt_exp.pdf', width=3.5, height=3)

par(mar=c(2.25,2.5,1.5,1), cex.axis=0.85)

plot(ev.fb$x, ev.fb$res, type='l', lwd=2, col='red',
     xlim=c(0, max(ev.fb$x)), ylim=c(0, 1),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, xaxs='i', yaxs='i')
polygon(c(ev.fb$x, rev(ev.fb$x)), c(ev.fb$res_10, rev(ev.fb$res_90)), 
        col=rgb(1,0,0,0.1), border=NA)
polygon(c(ev.fb$x, rev(ev.fb$x)), c(ev.fb$res_20, rev(ev.fb$res_80)), 
        col=rgb(1,0,0,0.2), border=NA)

legend('topleft', c('Median  | Residual RT |', 'Confidence Interval, 10%-90%\n(Smoothed)',
         'Confidence Interval, 20%-80%\n(Smoothed)'), 
       lty=c(1, NA, NA), lwd=c(2, NA, NA), pch=c(NA, 15, 15),
       col=c('red', rgb(1,0,0,0.1), rgb(1,0,0,0.3)),
       pt.cex=c(NA, 3, 3),
       bty='n', cex=0.75, y.intersp=1.5, inset=c(0.02, -0.06))

axis(1, at=seq(0, 220, by=40), tck=-0.01, mgp=c(0, 0.05, 0))
axis(2, at=seq(0, 1, by=0.2), tck=-0.01, mgp=c(0, 0.3, 0), las=1)

mtext('Experiment #, Sorted', side=1, cex=1, line=1.25)
mtext('| Residual RT |  (min)', side=2, cex=1, line=1.5)
mtext('Residual RT Varies by Experiment', side=3, cex=1, line=0.2, font=2)

dev.off()

## --------


ggsave(ggplot(ev.fb) +
  geom_ribbon(aes(x, ymin=res_10, ymax=res_90, 
                  fill="Confidence Interval, 10%-90%\n(Smoothed)")) +
  geom_ribbon(aes(x, ymin=res_20, ymax=res_80, 
                  fill="Confidence Interval, 20%-80%\n(Smoothed)")) +
  geom_path(aes(x, res, color="Mean |Residual RT|"), size=1) +
  scale_x_continuous(expand=c(0,0), limits=c(0, length(unique(ev.fb$`Raw file`))), 
                     breaks=seq(0, 280, by=40)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_manual(values=c(rgb(1,0,0,0.1), rgb(1,0,0,0.2))) +
  scale_color_manual(values=c('red'), labels=c("Median  | Residual RT |")) +
  guides(fill=guide_legend(order=2),
         color=guide_legend(order=1)) +
  labs(x="Experiment", y="| Residual RT | (min)",
       title='Residual RT Varies by Experiment',
       color=NULL, fill=NULL) +
  theme_bert() + theme(
    plot.margin=margin(0.2,0.1,0,0.1,'cm'),
    title=element_text(size=10),
    axis.text=element_text(size=10),
    plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'), hjust=0.5, face='bold'),
    legend.position=c(0.385, 0.9),
    legend.spacing.y=unit(0,'cm'),
    legend.margin=margin(0,0,0,0,'cm'),
    legend.key.height=unit(0.5,'cm'),
    legend.key.width=unit(0.5,'cm'),
    legend.text=element_text(size=10, margin=margin(0,0,1,0,'cm'))
  ), file='manuscript/Figs/fig_residual_rt_exp.pdf', width=3.5, height=3)
