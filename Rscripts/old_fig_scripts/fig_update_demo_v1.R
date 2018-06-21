## Update Demo - BAD ------

x <- seq(0, 250, by=0.1)
y1 <- dlnorm(x, meanlog=4, sdlog=0.5)
y2 <- dnorm(x, mean=120, sd=2)

# color palette
pal <- brewer.pal(4, 'Set1')

#plot(x, y2, type='l', col='blue')
#lines(x, y1, type='l', col='red')

density.plot <- 
  ggplot(data.frame()) +
  geom_path(aes(x=x, y=y2, color='b'), size=1) +
  geom_path(aes(x=x, y=y1, color='a'), size=1) +
  geom_segment(aes(x=90, y=0, xend=90, yend=0.25, color='c'), linetype='dashed', size=0.75) +
  geom_segment(aes(x=119, y=0, xend=119, yend=0.25, color='c'), linetype='dashed', size=0.75) +
  #annotate(geom='text', x=90, y=-0.05, label='Observation 1', size=5) +
  scale_color_manual(values=c(pal[1], pal[2], 'black'),
                     labels=c('Null RT Density (All RTs)', 
                              'Predicted RT Density\nfor Top Match', 
                              'Observed RT')) +
  scale_x_continuous(limits=c(0, 160), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.25), expand=c(0.01,0)) +
  labs(x=NULL, y='Density', color='Distribution') +
  theme_bert() +
  theme(
    plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'),
    #legend.position = c(0.1, 0.5),
    legend.position = c(0.23, 0.55),
    #legend.background = element_rect(fill='white', color=NULL, size=0.25),
    legend.background = element_blank(),
    legend.key.height = unit(1.1, 'cm'),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    #axis.line = element_line(size=0.25, color='black')
  )

den.plot <- ggplotGrob(density.plot)


bdf <- data.frame(
  a=factor(c('Spectra Only', 'Updated'), levels=c('Spectra Only', 'Updated')),
  b=as.numeric(c(1.9, 3))
)

bar.plot.1 <- ggplot(bdf) +
  #ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y='Error Prob. (PEP)\nfor Top Match') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.5)) +
  #annotate(geom='text', size=5, x=0.9, y=2.5, label='PSM B') +
  labs(title='PSM A') +
  theme_bert() %+replace% theme(
    plot.margin = unit(c(0.5,0.1,0.4,0.25), 'cm'),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    plot.title=element_text(size=16)
  )
bdf$b <- c(1.1, 0.5)
bar.plot.2 <- 
  ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y=' ') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.2)) +
  labs(title='PSM B') +
  theme_bert() %+replace% theme(
    plot.margin = unit(c(0.5,0.25,0.4,0.1), 'cm'),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    plot.title=element_text(size=16)
  )

bar.plot.1 <- ggplotGrob(bar.plot.1)
bar.plot.2 <- ggplotGrob(bar.plot.2)

#grid.arrange(den.plot.1, den.plot.2, bar.plot.1, bar.plot.2, 
#             ncol=2, nrow=2, widths=c(1,1), heights=c(3,2))

lay <- rbind(c(1,1,1,1,1,1),
             c(2,2,2,3,4,NA),
             c(5,5,5,6,6,6))
gs <- list(den.plot, 
           textGrob('Retention Time (mins)', hjust=0.6), 
           #textGrob('PSM A'), textGrob('PSM B'), 
           textGrob(''), textGrob(''), 
           bar.plot.1, bar.plot.2)
grid.arrange(grobs=gs, layout_matrix=lay, 
             widths=c(1,1,1,1,1,1), heights=c(1,0.1,1))

