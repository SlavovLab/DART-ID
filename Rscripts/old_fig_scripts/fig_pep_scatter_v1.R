## Fig 1A - PEP vs. PEP.new Scatter/Density -----

pdf(file='manuscript/Figs/Fig_1A.pdf', width=7, height=7, onefile=F)

p <- ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
  #ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
  stat_bin2d(bins=100, drop=TRUE, geom='tile', aes(fill=..count..)) +
  #stat_density2d(aes(alpha=..level.., fill=..level..), bins=9, geom='polygon') +
  geom_abline(slope=1, intercept=0, color='red') +
  scale_fill_gradientn(colors=bHeatmap, 
                       values=c(0, 0.05, 0.1, 0.2, 0.5, 1)
                       #labels=c(0, 500, 1000, 1500, 2000)) +
  )+
  scale_alpha_continuous(guide=FALSE) +
  scale_x_log10(limits=c(1e-10, 1), 
                breaks=logseq(1e-10, 1, 6), labels=fancy_scientific) +
  #breaks=logseq(1e-10, 1, 6)) +
  scale_y_log10(limits=c(1e-10, 1), 
                breaks=logseq(1e-10, 1, 6), labels=fancy_scientific) +
  #breaks=logseq(1e-10, 1, 6)) +
  theme_bert() +
  theme(legend.position=c(0.9, 0.3),
        legend.key.height = unit(0.05, 'npc'),
        legend.text = element_text(colour = "black", size = 12),
        plot.margin=unit(c(0,0,0.5,0.5), 'cm')) +
  labs(fill='Count', x='Spectral PEP', y='Updated PEP')

density.top <- ggplot(ev.f, aes(PEP)) + 
  #ggplot(ev.f, aes(PEP)) + 
  geom_histogram(bins=30, aes(y=..density..), color='black', fill='white') +
  stat_density(adjust=8, geom='path', color='red') + 
  scale_x_log10(limits=c(1e-10, 1)) +
  labs(y='Density') +
  theme_bert() +
  theme(axis.text=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x=element_line(color='ivory3', size=0.5),
        panel.grid.minor.x=element_line(color='ivory2', size=0.5),
        plot.margin=unit(c(0.5,0.5,0,0.5), 'cm'))
density.right <- ggplot(ev.f, aes(PEP.new)) + 
  geom_histogram(bins=30, aes(y=..density..), color='black', fill='white') +
  stat_density(adjust=5, geom='path', color='red') + 
  scale_x_log10(limits=c(1e-10, 1)) + 
  labs(y='Density') +
  coord_flip() +
  theme_bert() +
  theme(axis.text=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_line(color='ivory3', size=0.5),
        panel.grid.minor.y=element_line(color='ivory2', size=0.5),
        plot.margin=unit(c(0.5,0.5,0.5,0), 'cm'))

bTitle <- textGrob(label='[Method Name] PEP Shift', 
                   x=unit(1, 'cm'),
                   just=c('left', 'centre'),
                   gp=gpar(fontsize=24))

main.grob <- ggplotGrob(p)
den.top.grob <- ggplotGrob(density.top)
den.right.grob <- ggplotGrob(density.right)

den.top.grob$widths <- main.grob$widths
den.right.grob$heights <- main.grob$heights

grid.arrange(#bTitle,        nullGrob(),
  den.top.grob,  nullGrob(),
  main.grob,     den.right.grob, 
  #nrow=3, ncol=2, widths=c(5, 1), heights=c(1, 2, 7))
  nrow=2, ncol=2, widths=c(5, 1), heights=c(1, 5))

dev.off()
