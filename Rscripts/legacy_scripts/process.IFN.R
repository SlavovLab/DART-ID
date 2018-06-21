# build IFN data set
library(tidyverse)

prefix <- '/Volumes/GoogleDrive/My Drive/Ribosomes_MS_Data/'
suffix <- 'evidence.txt'

files = c(
  'IFN1/IFN1A/txt/',
  'IFN1/IFN1B/txt/',
  'IFN1/IFN1C/txt/',
  'IFN1/IFN1D/txt/',
  'IFN1/IFN1E/txt/',
  'IFN2/IFN2A/txt/',
  'IFN4/IFN4A/txt/',
  'IFN5/IFN5ABCtxt/',
  'IFN7/IFN7A/txt/',
  'IFN7/IFN7BCD/txt/'
)
files <- paste0(prefix, files, suffix)

ev <- data.frame()
for(file in files) {
  ev.a <- read_tsv(file)
  ev <- rbind(ev, ev.a[,c(
    'Sequence', 'Modified sequence', 'Proteins', 'Leading razor protein', 'Type',
    'Raw file', 'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity'
    # ,'Reporter intensity corrected 0', 'Reporter intensity corrected 1',
    # 'Reporter intensity corrected 2', 'Reporter intensity corrected 3',
    # 'Reporter intensity corrected 4', 'Reporter intensity corrected 5',
    # 'Reporter intensity corrected 6', 'Reporter intensity corrected 7',
    # 'Reporter intensity corrected 8', 'Reporter intensity corrected 9'
  )])
}

ev <- ev %>%
  mutate(
    `Peptide ID`=as.numeric(as.factor(ev$`Modified sequence`)),
    `Sequence ID`=as.numeric(as.factor(ev$Sequence))
  ) %>%
  arrange(`Peptide ID`) %>%
  mutate(id=seq(1,nrow(ev)),
         `Mod. peptide ID`=`Peptide ID`)

cor.mat <- filter.exps(ev, pep.thresh=1e-2)
plot.cor.mat(cor.mat, show.text=T)

# remove 14
ev <- ev[!ev$`Raw file`==raw.files[14],]

write_tsv(ev, 'dat/ev.IFN.txt')


# run alignment
ev <- read_tsv('dat/ev.IFN.txt')

align.rt(ev, 
         pep_thresh=0.1,
         prior.iters=20,
         pars.out='dat/params.IFN.RData',
         figs.out='align_IFN_PEP_0_1')
