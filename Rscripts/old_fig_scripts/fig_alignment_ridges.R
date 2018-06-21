## -----

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'PEP.new')) %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=log10(PEP/PEP.new)) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))


## Alignment Ridges -----

ev.f <- ev %>%
  select(c('Sequence', 'Leading razor protein', 'Raw file', 'PEP', 'PEP.new', 'Peptide ID', 
           'Retention time', 'muijs', 'sigmas')) %>%
  filter(!is.na(PEP.new)) %>%
  #filter(PEP < 0.05) %>%
  filter(!grepl('REV*', `Leading razor protein`)) %>% # Remove Reverse matches
  filter(!grepl('CON*',`Leading razor protein`)) # Remove Contaminants

# clean up raw file names and extract an experiment ID from it
# need the experiment ID (19A, 30B, etc.), to match with sample metadata from the
# experiment description excel sheet/.csv
ev.f$file <- clean.file.name(ev.f$`Raw file`)
ev.f$exp <- str_extract(ev.f$file, '[1-9][0-9][A-Z](\\d)?')

psm.count <- ev.f %>%
  group_by(`Peptide ID`) %>%
  summarise(Sequence=unique(Sequence), 
            Protein=unique(`Leading razor protein`), 
            count=length(`Peptide ID`)) %>%
  arrange(desc(count))



## ------

library(ggridges)
library(rmutil)


i <- 5
pep.id <- psm.count$`Peptide ID`[i]
ev.p <- ev.f %>% 
  filter(`Peptide ID`==pep.id)

# take the top [num.files]
num.files <- 5

alignment.dat <- ev.p %>%
  group_by(exp) %>%
  summarise(count=length(exp), 
            mu=first(unique(muijs)), 
            sigma=first(unique(sigmas))) %>%
  arrange(desc(count)) %>%
  top_n(num.files, wt=count)

files <- alignment.dat$exp[1:num.files]

density.dat <- data.frame(x=numeric(), y=numeric(), exp=character())
range <- seq(0,250,by=0.5)
for(j in 1:num.files) {
  density.dat <- rbind(density.dat, data.frame(
    x=as.numeric(range),
    y=as.numeric(dlaplace(range, m=alignment.dat$mu[j], s=alignment.dat$sigma[j])),
    exp=as.character(files[j])
  ))
}

ev.p %>% 
  filter(exp %in% files) %>%
  filter(PEP.new < 0.05) %>%
  ggplot() +
  geom_density_ridges(aes(x=`Retention time`, y=exp, group=exp), 
                      rel_min_height=0.01) +
  geom_path()
#geom_hline(yintercept=2, color='red') +
#scale_x_continuous(limits=c(40, 150)) +
labs(title=unique(ev.p$Sequence))

