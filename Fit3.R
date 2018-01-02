library(tidyverse)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3.txt')

ev <- ev %>%
  mutate(eRT=`Retention time`-muijs)

## look at PSMs which get severely downgraded by the fit
downgrades <- ev %>%
  filter(PEP < 1e-5) %>%
  filter(PEP.new > 1e-1) %>%s
  mutate(dRT=`Retention time`-muijs) %>%
  pull(`Peptide ID`)

sum(ev$`Peptide ID` %in% downgrades)

for(pid in downgrades) {
  print(pid)
  ev.d <- ev %>% filter(`Peptide ID`==pid)
  #cat(id, nrow(ev.d), '\n')
}

