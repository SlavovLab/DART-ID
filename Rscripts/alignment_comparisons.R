library(tidyverse)
#library(gdata)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ev_updated.txt')

ev.f <- ev %>%
  filter(!is.na(pep_new)) %>%
  filter(!grepl('SQC9', `Raw file`)) %>%
  filter(PEP < 0.01)

## create ssrcalc input -----

length(unique(ev.f %>% filter(PEP < 0.01) %>% pull(Sequence)))
write(unique(ev.f %>% filter(PEP < 0.01) %>% pull(Sequence))[1:3000], '/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrcalc_sequence_list_1.txt')
write(unique(ev.f %>% filter(PEP < 0.01) %>% pull(Sequence))[3001:length(unique(ev.f %>% filter(PEP < 0.01) %>% pull(Sequence)))], '/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrcalc_sequence_list_2.txt')

# http://hs2.proteome.ca/SSRCalc/SSRCalcQ.html
# separation system: 100Å C18 column, 0.1% Formic Acid 2015
# label deltas: TMT
# cys protection: Free Cysteine	

# no retention times, or HI inputted. no chart produced.

## load ssrcalc output ------

hi_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA323_313_list1.txt')
hi_dat <- rbind(hi_dat, read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA7794_648_list2.txt'))

# append predicted HI to DART dataframe
ev.f$HI <- hi_dat$`HI (pred)`[match(ev.f$Sequence, hi_dat$Sequence)]

ev.f %>% filter(PEP < 0.01) %>% #filter(`Raw file`=='180424S_X_IFN6J') %>%
ggplot(aes(x=HI, y=`Retention time`)) +
  geom_point(alpha=0.5)

# for each experiment, build linear relationship between HI and RT,
# then get errors for the line and the points
# exps <- sort(unique(ev.f$`Raw file`))
# ssrcalc_error <- sapply(exps, function(exp){
#   ev.a <- ev.f %>% filter(PEP < 0.01) %>% filter(`Raw file`==exp) %>% filter(!is.na(HI))
#   coefs <- lm(ev.a$`Retention time` ~ ev.a$HI)$coefficients
#   ev.a$`Retention time` - ((ev.a$HI * coefs[2]) + coefs[1])
# })

# are the means and sds roughly the same?
# sapply(ssrcalc_error, function(err) {
#   sd(err)
# })

# global error
ssrcalc_coefs <- lm(ev.f$`Retention time` ~ ev.f$HI)$coefficients
ssrcalc_error <- ev.f$`Retention time` - ((ev.f$HI * ssrcalc_coefs[2]) + ssrcalc_coefs[1])

## load biolccc input ---------

View(ev.f %>%
  group_by(`Raw file`) %>%
  summarise(n=max(`Retention time`)))

seqs <- sort(unique(ev.f %>% 
  filter(PEP < 0.01) %>% 
  filter(!grepl('95', `Raw file`)) %>%
  filter(`Retention time` < 48) %>%
  pull(`Modified sequence`)))

seqs <- gsub('([A-Z])_', '\\1', seqs, perl=T)
seqs <- gsub('_\\(ac\\)', 'Ac-', seqs)
seqs <- gsub('\\(ox\\)', '', seqs)
seqs <- gsub('_([A-Z])', '\\1', seqs, perl=T)

write(seqs, '/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/biolccc_sequence_list.txt')

# Column length, mm	250.0
# Column internal diameter, mm	0.075
# Packing material pore size, A	130.0
# Initial concentration of component B, %	5.0
# Final concentration of component B, %	35.0
# Gradient time, min	48.0
# Delay time, min	0.0
# Flow rate, ml/min	0.0001
# ACN concentration in component A, %	0.0
# ACN concentration in component B, %	80.0
# Solid/mobile phase combination	RP/ACN+FA
# no cysteine carboxyaminomethylation
# http://www.theorchromo.ru/

## biolccc output --------

biolccc_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/biolccc_output.txt')

# i know acetylated and oxidized sequences will be thrown out, but don't need them for a quick comparison of the method
ev.f$biolccc_rt <- biolccc_dat$`Retention time, min`[match(ev.f$Sequence, biolccc_dat$Sequence)]

# for each experiment, build linear relationship between predicted RT and RT,
# then get errors for the line and the points
# exps <- sort(unique(ev.f$`Raw file`))
# biolccc_error <- sapply(exps, function(exp){
#   ev.a <- ev.f %>% 
#     filter(PEP < 0.01) %>% 
#     filter(`Raw file`==exp) %>% 
#     filter(!is.na(biolccc_rt))
#   coefs <- lm(ev.a$`Retention time` ~ ev.a$biolccc_rt)$coefficients
#   ev.a$`Retention time` - ((ev.a$biolccc_rt * coefs[2]) + coefs[1])
# })
# 
# sapply(biolccc_error, function(err) {
#   sd(err)
# })

ev.a <- ev.f %>% filter(!is.na(biolccc_rt))
biolccc_coefs <- lm(ev.a$`Retention time` ~ ev.a$biolccc_rt)$coefficients
biolccc_error <- ev.f$`Retention time` - ((ev.f$biolccc_rt * biolccc_coefs[2]) + biolccc_coefs[1])

## elude input -------

ev.f <- ev %>% 
  filter(!grepl('95', `Raw file`)) %>% 
  filter(!is.na(pep_new)) %>%
  filter(PEP < 0.01) %>% 
  filter(Modifications=='Unmodified')

# split dataset in half for training and test sets.
seqs <- sort(unique(ev.f$Sequence))
set.seed(1)
training_seqs <- sort(sample(seqs, floor(length(seqs) / 2)))
testing_seqs <- seqs[!seqs %in% training_seqs]

ev.training <- ev.f %>% filter(Sequence %in% training_seqs) %>% select(c('Sequence', 'Retention time'))
write_tsv(ev.training, path='/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_training.txt', col_names=F)

ev.testing <- ev.f %>% filter(Sequence %in% testing_seqs) %>% select(c('Sequence', 'Retention time'))
write_tsv(ev.testing, path='/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_testing.txt', col_names=F)

# /usr/local/bin/elude.app/Contents/MacOS/elude -t /gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_training.txt -e /gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_testing.txt -s /gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude.model -o /gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_out.txt -y -g -v 10

## elude output ---------

elude_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_out.txt', skip=2)
elude_dat <- elude_dat %>% 
  arrange(Peptide) %>% group_by(Peptide) %>%
  summarise(RT=unique(Predicted_RT))

ev.f$elude_RT <- elude_dat$RT[match(ev.f$Sequence, elude_dat$Peptide)]

# for each experiment, build linear relationship between predicted RT and RT,
# then get errors for the line and the points
# exps <- sort(unique(ev.f$`Raw file`))
# elude_error <- sapply(exps, function(exp){
#   ev.a <- ev.f %>% filter(PEP < 0.01) %>% filter(`Raw file`==exp) %>% filter(!is.na(elude_RT))
#   if(nrow(ev.a) == 0) return(c(0))
#   coefs <- lm(ev.a$`Retention time` ~ ev.a$elude_RT)$coefficients
#   ev.a$`Retention time` - ((ev.a$elude_RT * coefs[2]) + coefs[1])
# })
# 
# sapply(elude_error, function(err) {
#   sd(err)
# })

ev.a <- ev.f %>% filter(!is.na(elude_RT))
elude_coefs <- lm(ev.a$`Retention time` ~ ev.a$elude_RT)$coefficients
elude_error <- ev.f$`Retention time` - ((ev.f$elude_RT * elude_coefs[2]) + elude_coefs[1])


# MaxQuant match-between-runs ---------------------------------------------

ev.f %>% #filter(exp_id %in% c(5, 7)) %>%
  ggplot(aes(`Calibrated retention time`, `Retention time`)) +
  geom_point(alpha=0.5) +
  labs(x='Calibrated retention time (mins)', y='Retention time (min)') +
  theme(axis.text=element_text(size=12), title=element_text(size=18))

exps <- sort(unique(ev.f$`Raw file`))
ev.mq <- data.frame()
for(exp in exps) {
  ev.a <- ev.f %>% filter(PEP < 0.01) %>% filter(`Raw file`==exp) %>%
    filter(!is.na(`Calibrated retention time`))
  if(nrow(ev.a) == 0) return(c(0))
  coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
  #ev.a$`Retention time` - ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.mq <- rbind(ev.mq, ev.a)
}
#sapply(mqmbt_error, function(err) {
#  sd(err)
#}
plot(ev.mq$`Retention time`, ev.mq$RT_calibrated)
hist(ev.mq$`Retention time`-ev.mq$RT_calibrated)

ev.a <- ev.f %>% filter(PEP < 0.01)
coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
plot(ev.a$`Retention time`[1:1e4], ev.a$RT_calibrated[1:1e4])
abline(a=0, b=1, col='red')

#ev.a <- ev.f %>% filter(!is.na(`Calibrated retention time`))
#mqmbr_coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
#mqmbr_error <- ev.f$`Retention time` - ((ev.f$`Calibrated retention time` * mqmbr_coefs[2]) + mqmbr_coefs[1])
mqmbr_error <- ev.mq$`Retention time`-ev.mq$RT_calibrated

## DART ---------

# exps <- sort(unique(ev.f$`Raw file`))
# dart_error <- sapply(exps, function(exp) {
#   ev.a <- ev.f %>% filter(PEP < 0.01) %>% filter(`Raw file`==exp)
#   ev.a$`Retention time` - ev.a$muij
# })
# 
# sapply(dart_error, function(err) {
#   #sd(err)
#   mean(err)
# })

dart_error <- ev.f$`Retention time` - ev.f$muij


# iRT ---------------------------------------------------------------------

# fake data
# ev.f$irt_rt <- ev.f$`Retention time` + rnorm(nrow(ev.f), mean=0, sd=3)
# irt_error <- ev.f$`Retention time` - ev.f$irt_rt

#df_irt <- gdata::read.xls('dat/irt.xls', verbose=T, method='tab')
df_irt <- read_tsv('dat/irt.txt')

df_irt_f <- df_irt %>% 
  filter(!is.na(PP.DeltaRT)) %>%
  filter(PEP.StrippedSequence %in% ev.f$Sequence) %>%
  arrange(PEP.StrippedSequence) %>%
  group_by(PEP.StrippedSequence, PP.DeltaRT) %>%
  summarise(RT=unique(PP.EmpiricalRT),
            predRT=unique(PP.RTPredicted)) %>%
  rename(Sequence=PEP.StrippedSequence, dRT=PP.DeltaRT)


# combine errors into DF ---------------------------------------------------------

# order: ssrcalc, elude, biolccc, irt, mqmbt, dart
types <- c('SSRCalc', 'ELUDE', 'BioLCCC', 'MaxQuant MBR', 'DART-ID')
error_df <- data.frame(
  pred=as.numeric(c(ev.f$HI, ev.f$elude_RT, ev.f$biolccc_rt, 
                    #ev.f$irt_rt, 
                    ev.mq$`Calibrated retention time`, ev.f$muij)),
  pred_adjusted=as.numeric(c((ev.f$HI*ssrcalc_coefs[2]) + ssrcalc_coefs[1],
                             (ev.f$elude_RT*elude_coefs[2]) + elude_coefs[1],
                             (ev.f$biolccc_rt*biolccc_coefs[2]) + biolccc_coefs[1],
                             #ev.f$irt_rt,
                             ev.mq$RT_calibrated,
                             ev.f$muij)),
  obs=as.numeric(c(ev.f$`Retention time`, ev.f$`Retention time`, ev.f$`Retention time`,
                   ev.mq$`Retention time`, ev.f$`Retention time`)),
  error=as.numeric(c(ssrcalc_error, elude_error, biolccc_error, 
                     #irt_error,
                     mqmbr_error, dart_error)),
  type=factor(rep(types, each=nrow(ev.f)), levels=types)
)

# add iRT data

error_df <- rbind(error_df, data.frame(
  pred=as.numeric(df_irt_f$predRT),
  pred_adjusted=as.numeric(df_irt_f$predRT),
  obs=as.numeric(df_irt_f$RT),
  error=as.numeric(df_irt_f$dRT),
  type='iRT'
))

error_df$type <- factor(error_df$type, levels=types <- c('SSRCalc', 'ELUDE', 'BioLCCC', 'iRT', 'MaxQuant MBR', 'DART-ID'))


# save data ---------------------------------------------------------------

save(error_df, file='dat/error_df_20180716.RData')


# summary stats -----------------------------------------------------------

error_df %>%
  group_by(type) %>%
  summarise(sd=sd(error, na.rm=T),
            mean=mean(abs(error), na.rm=T),
            median=median(abs(error), na.rm=T))
