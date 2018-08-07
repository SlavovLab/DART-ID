library(tidyverse)
library(ggridges)
source('Rscripts/lib.R')

ev <- read_tsv('~/git/RTLib/Alignments/PS_60min_20180627_2/ev_updated.txt')

## ----------

ev <- ev %>%
  mutate(residual=`Retention time`-muij) %>%
  mutate(residual_abs=abs(residual))

ggplot(ev) +
  #geom_boxplot(aes(group=`Raw file`)) +
  #geom_violin(aes(group=`Raw file`)) +
  geom_density_ridges(aes(x=residual, y=`Raw file`, group=`Raw file`),
                      rel_min_height=0.01, scale=3) +
  scale_x_continuous(limits=c(-0.5, 0.5)) +
  labs(x='Residual RT (min)', title='PS - 60 min runs') +
  theme_ridges() +
  theme(
    axis.text.y=element_text(size=8)
  )

mean(ev$residual_abs, na.rm=T) * 60

## -----------

ev %>%
  group_by(`Raw file`) %>%
  summarise(m=mean(residual_abs, na.rm=T)) %>%
  pull(m)

## --------
ev.f <- ev %>%
  filter(grepl('(ph)', `Modified sequence`))

sum(ev.f$`Spectra PEP` < 0.01)
sum(ev.f$PEP < 0.01)

