## mixture model demo ------

x <- seq(0,60,by=0.1)
#y1 <- dlnorm(x, meanlog=4.663, sdlog=0.5089)
y1 <- dnorm(x, mean=38, sd=17)
y2 <- dnorm(x, mean=20, sd=1.78)
#y3 <- dnorm(x, mean=80, sd=2.3)
#y4 <- dnorm(x, mean=120, sd=2)

#plot(x, y2, 'l', col='red')
#lines(x,y1,'l', col='black')

#p <- ggplot(data.frame(x,y1,y2,y3,y4)) +
p <- ggplot(data.frame(x,y1,y2)) +
  geom_area(aes(x=x, y=y1), fill='red', alpha=0.3) + 
  geom_path(aes(x=x, y=y1), color='red', size=0.4) +
  geom_area(aes(x=x, y=y2), fill='blue', alpha=0.3) + 
  geom_path(aes(x=x, y=y2), color='blue', size=0.4) +
  #geom_path(aes(x=x, y=y3), color='blue', size=0.4) +
  #geom_path(aes(x=x, y=y4), color='blue', size=0.4) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.24), expand=c(0,0)) +
  labs(x='Retention Time', y='Density') +
  theme_bw() %+replace% theme(
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    axis.title.x = element_text(family='Helvetica', size=6),
    axis.title.y = element_text(family='Helvetica', size=6, angle=90),
    panel.grid.minor=element_blank()
    #panel.grid=element_blank()
  )
ggsave('manuscript/Figs/mixture_model_demo.pdf', plot=p, 'pdf', width=4, height=2.5, units='cm')


