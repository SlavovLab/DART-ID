## PEP Fold Change Line Plot -------

source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.txt')
ev.perc <- read_tsv('dat/ev.perc_181217_noPEP.txt')
ev.perc.nodoc <- read_tsv('dat/ev.perc_181217_nodoc.txt')

## -----

exps <- list(RTLib=ev,
             Percolator=ev.perc)
#Percolator_NoDoc=ev.perc.nodoc)

df <- fold.change.comp(exps, num.steps=100)
df <- fold.change.comp(exps, begin=0, end=1, num.steps=30, log=F)

ggplot(df, aes(x=x, y=PEP, color=Method)) +
  geom_path() +
  scale_x_log10() +
  annotation_logticks(sides='b') +
  #scale_y_continuous(limits=c(0.95,2.25), breaks=c(1, 1.25, 1.5, 1.75, 2, 2.25)) +
  labs(x='PEP Threshold', y='Fold Change Increase in IDs',
       title=paste0('Fold Change Increase of PSM IDs\n',
                    '= #Adjusted PEPs / #Original PEPs above PEP Threshold'))

pdf('manuscript/Figs/Fig_1B.pdf')

ggplot(df, aes(x=x, y=PEP, color=Method)) +
  #geom_hline(yintercept=1, color='red', linetype='longdash') +
  geom_path() +
  #geom_point() +
  scale_x_log10(name='PEP Interval', breaks=logseq(1e-5, 1, 6)) +
  #labels=fancy_scientific) +
  annotation_logticks(sides='b') +
  #scale_x_continuous(breaks=seq(0, 1, by=0.2)) +
  scale_y_continuous(breaks=seq(0.4,2,by=0.2)) +
  labs(x='PEP Threshold', 
       y=parse(text='frac(\'# Updated PEP < Threshold\', \'# Spectral PEP < Threshold\')'),
       #y='asdf',
       #title=parse(text="frac(\'#Adjusted PEPs\', \'#Original PEPs above PEP Threshold\')")) +
       title=paste0("Fold Change Increase of Confident PSMs")) +
  theme_bert()

dev.off()
