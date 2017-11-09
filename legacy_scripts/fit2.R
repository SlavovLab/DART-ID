library(pracma)
library(tidyverse)
source('parse.ev.adj.R')
source('adjust.pep.ali.R')

ev.pep <- parse.ev.adj('dat/ev.adj.elite.txt')
ev.fit2 <- parse.ev.adj('dat/ev.adj.Fit2.txt')


## compare fit2 to normal fit
ev.pep$Method <- 'Peptide-Centric (STAN)'
ev.fit2$Method <- 'Fit 2 (STAN w/ updated sigmas)'
df <- rbind(
  ev.pep[,c('dPEP', 'dRT', 'Retention time', 'RT.new', 'Method', 'PEP', 'PEP.new')],
  ev.fit2[,c('dPEP', 'dRT', 'Retention time', 'RT.new', 'Method', 'PEP', 'PEP.new')]
)
df <- df[sample.int(nrow(df), size=1e5),]
df <- df %>%
  mutate(method.id=Method) %>%
  mutate_at('method.id', funs(as.numeric(as.factor(.)))) %>%
  mutate(PEP.class.new=ifelse(PEP.new < 5e-2, 'PEP.new < 0.05', 'PEP.new > 0.05')) %>%
  mutate(PEP.class=ifelse(PEP < 5e-2, 'PEP < 0.05', 'PEP > 0.05')) %>%
  mutate(PEP.new = ifelse(PEP.new < 1e-4, 1e-4, PEP.new))
df %>%
  filter(!is.na(PEP.new)) %>%
  #filter(method.id==1) %>%
  #ggplot(aes(x=`Retention time`, y=dPEP)) +
  #ggplot(aes(x=`Retention time`, y=(dRT))) +
  ggplot(aes(x=`Retention time`, y=RT.new)) +
  #geom_point(aes(color=PEP.new), alpha=0.1) +
  geom_point(aes(color=PEP), alpha=0.1) +
  #geom_point(aes(color=abs(dPEP)), alpha=0.1) + 
  #scale_color_gradient(low='red', high='green', trans='log', 
                       #limit=c(1e-10, 1), breaks=c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1)) +
  scale_color_gradientn(trans='log', limit=c(1e-4, 1), na.value='#006837',
                        breaks=c(1e-4, 1e-2, 1e-1, 1),
                        #colors=brewer.pal(name='RdYlGn', n=11), 
                        colors=c('#006837', '#1A9850', '#D9EFB8', '#A50026'),
                        values=c(1e-4, 1e-2, 1e-1, 1)) +
  #stat_density2d(aes(fill=..level..), geom='polygon', n=100) +
  #geom_hline(yintercept=0, color='red', linetype='longdash') +
  geom_abline(intercept=0, slope=1, color='black', linetype='longdash') +
  #facet_grid(~Method) +
  #facet_grid(PEP.class.new ~ Method) +
  facet_grid(PEP.class ~ Method) +
  labs(
       #fill='Density', 
       x='Retention Time (mins)', 
       #y='dPEP (PEP.new - PEP)',
       #y='Fit Residuals',
       y='RT.new (mins)',
       title='RTLib: RT vs. RT.new')
  #geom_point(alpha=0.01) +
#scale_y_log10()
#scale_y_continuous()

## look at it experiment by experiment


