source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3.txt')
ev <- read_tsv('dat/ev.adj.Fit3.RTL.txt')

## alignment demo -----

x <- seq(0,100,by=0.1)
split_point <- 40
beta_0 <- 10
beta_1 <- 0.8
beta_2 <- 2
y <- vector(length=length(x), mode='numeric')
y[x <= split_point] <- beta_0 + (x[x <= split_point] * beta_1)
y[x > split_point] <- beta_0 + (split_point * beta_1) + ((x[x > split_point]-split_point) * beta_2)

#plot(x, y, 'l', col='blue', lwd=2)
#abline(v=split_point, col='red', lwd=2, lty='dashed')

p <- ggplot(data.frame(x,y), aes(x,y)) +
  geom_path(color='blue', size=0.5) +
  geom_vline(xintercept=split_point, color='red', size=0.5) +
  scale_x_continuous(expand=c(0,0)) +
  #labs(x=parse(text='mu[i]'), y=parse(text='rho[i][j][k]')) +
  labs(x=NULL, y=NULL) +
  theme_bw() %+replace% theme(
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    #panel.grid=element_blank()
  )
ggsave('manuscript/Figs/alignment_demo.pdf', plot=p, 'pdf', width=3, height=3, units='cm')
