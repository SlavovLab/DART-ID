## comparison to match between runs --------

#ev_mq <- read_tsv("/gd/MS/SCoPE/SQC/SQC87_95_match_between_runs/evidence.txt")
ev_mqi <- read_tsv("/gd/MS/SCoPE/SQC/SQC_67_95_Varied/evidence.txt")
#ev_dt <- read_tsv("Alignments/SQC87_95_20180611_2/ev_updated.txt")
ev_dti <- read_tsv("Alignments/SQC_varied_20180613_5/ev_updated.txt")

## --------

ev_mq <- ev_mqi %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev_mqi)))[order(order(PEP))]) %>%
  filter(qval < 0.01) %>%
  #filter(!Type=='MULTI-MATCH') %>%
  group_by(Sequence) %>%
  summarise(res_mean=mean(`Calibrated retention time`),
            res_sd=sd(`Calibrated retention time`),
            num=n()) %>%
  filter(!is.na(res_sd)) %>%
  filter(num > 10)

ev_dt <- ev_dti %>%
  #mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  #mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev_dti)))[order(order(PEP))]) %>%
  mutate_at('pep_updated', funs(ifelse(. > 1, 1, .))) %>%
  mutate(qval=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(ev_dti)))[order(order(pep_updated))],
         dPEP=log10(PEP/pep_updated)) %>%
  filter(!is.na(muij)) %>%
  filter(qval < 0.01) %>%
  group_by(Sequence) %>%
  summarise(res_mean=mean(`Retention time`-muij),
            res_sd=sd(`Retention time`-muij),
            num=n()) %>%
  filter(!is.na(res_sd)) %>%
  filter(num > 10)

peps <- intersect(ev_mq$Sequence, ev_dt$Sequence)

ev_mq <- ev_mq %>% filter(Sequence %in% peps)
ev_dt <- ev_dt %>% filter(Sequence %in% peps)

## --------

ev_mqi %>% 
  mutate(exp_id=as.numeric(as.factor(`Raw file`))) %>%
  filter(exp_id==3) %>%
  #filter(exp_id %in% c(1, 2, 3, 4)) %>%
  #sample_n(1e3) %>% 
ggplot(aes(x=`Retention time`, y=`Calibrated retention time`)) + 
  geom_point(alpha=0.2)

seq_freq <- as.data.frame(table(paste(ev_mqi$`Modified sequence`, ev_mqi$Charge))) %>% 
  arrange(desc(Freq))

ev_mqi %>%
  filter(paste(`Modified sequence`, Charge) == as.character(seq_freq[4,]$Var1)) %>%
ggplot(aes(x=`Retention time`, y=`Calibrated retention time`)) + 
  geom_point(alpha=0.2)


## --------

plot(density(ev_dt$res_sd, na.rm=T, adjust=2),
     col=av[2], lwd=2,
     xlim=c(0, 2), ylim=c(0, 4),
     xlab='Residual RT (min)', ylab='Density',
     main='asdf')
lines(density(ev_mq$res_sd, na.rm=T, adjust=2),
      col=av[1], lwd=2)

## --------

plot(ev_mq$res_sd, ev_dt$res_sd, pch=16, col=rgb(0,0,0,0.1),
     xlim=c(0, 5), ylim=c(0, 5),
     xlab='std RT - MaxQuant (min)', ylab='std RT - DART-ID (min)',
     main=paste0('RT Alignment Residuals\n','n=',nrow(ev_mq)))
abline(a=0, b=1, col='red', lwd=2)

## -------

ev_dt <- ev_dti %>%
  filter(!is.na(muij)) %>%
  filter(!Type=='MSMS') %>%
  filter(pep_updated < 0.01) %>%
  mutate(residual=`Retention time`-muij) %>%
  mutate(pep_id=paste(`Modified sequence`,Charge,`Raw file`,sep="-")) %>%
  arrange(Sequence)

ev_mq <- ev_mqi %>%
  #filter(PEP < 0.01) %>%
  filter(!is.na(`Match time difference`)) %>%
  mutate(pep_id=paste(`Modified sequence`,Charge,`Raw file`,sep="-"))

# get the intersection between these two sets
pep_exp_pairs = intersect(ev_dt$pep_id, ev_mq$pep_id)
ev_dt <- ev_dt %>% filter(pep_id %in% pep_exp_pairs)
ev_mq <- ev_mq %>% filter(pep_id %in% pep_exp_pairs)

## -------

plot(density(ev_dt$residual, n=1e4), col=av[2], lwd=2,
     xlab='Residual RT (min)', ylab='Density',
     main='RT Alignment Comparison',
     xlim=c(-1, 1))
lines(density(ev_mq$`Match time difference`, n=1e4), col=av[1], lwd=2)
legend('topleft', c("MaxQuant MBR", "DART-ID"),
       col=c(av[1], av[2]), lty=c(1, 1), lwd=c(2, 2),
       cex=0.8, bty='n', y.intersp=1.4)

## -------

plot(density(abs(ev_dt$residual), n=1e4), xlim=c(0, 1),ylim=c(0, 6))
lines(density(abs(ev_mq$`Match time difference`), n=1e4))

boxplot(abs(ev_mq$`Match time difference`), abs(ev_dt$residual), 
        ylim=c(0, 0.4), outline=F,
        names=c("MaxQuant MBR", "DART-ID"),
        col=c(av[1],av[2]))

## --------

plot(
  x=seq(-1, 1, by=1e-3),
  y=ecdf(ev_mq$`Match time difference`)(seq(-1, 1, by=1e-3)),
  type='l', col=av[1], lwd=2, 
  xlab='', ylab=''
  #xaxt='n', yaxt='n'
)
lines(
  x=seq(-1, 1, by=1e-3),
  y=ecdf(ev_dt$residual)(seq(-1, 1, by=1e-3)),
  col=av[2], lwd=2
)
title(x="Residual RT (min)", y="Density")
legend('topleft', c("MaxQuant MBR", "DART-ID"),
       col=c(av[1], av[2]), lty=c(1, 1), lwd=c(2, 2),
       cex=0.8, bty='n', y.intersp=1.4)


## -------------



