## dPEP Ridges ------

ev.f <- ev %>%
  filter(!is.na(PEP.new)) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'PEP.new', 'Retention time', 'muijs')) %>%
  mutate(PEP=ifelse(PEP > 1, 1, PEP), 
         PEP=ifelse(PEP == 0, .Machine$double.xmin, PEP)) %>%
  mutate(dPEP=log2(PEP/PEP.new),
         dRT=log10(abs(`Retention time` - muijs))) %>%
  mutate_at(vars(dPEP, dRT), funs(ifelse(is.infinite(.), NA, .))) %>%
  #mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))
  mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.1)))
#mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.2)))

## dPEP/dRT Ridges - Plotting ----

ggplot(ev.f) +
  #geom_density(aes(dPEP))
  geom_density_ridges(aes(x=dRT, y=bin, group=bin), rel_min_height=0.01, scale=1.2) +
  #stat='binline', bins=60) +
  geom_vline(xintercept=0, color='red') +
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, by=1)) +
  scale_y_discrete(expand=c(0.05,0)) +
  labs(x='dRT (min)', y='Spectral PEP') +
  theme_ridges() +
  theme(
    axis.title.x = element_text(size=16, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=16, hjust=0.5, vjust=0.5)
  )