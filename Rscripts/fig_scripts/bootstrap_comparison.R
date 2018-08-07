library(tidyverse)
source('Rscripts/lib.R')

ev.a <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180802_5exp_point/ev_updated.txt")
ev.b <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180802_5exp_non-parametric/ev_updated.txt")
ev.c <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180802_5exp_parametric/ev_updated.txt")
ev.d <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180803_5exp_parametric_mixture_v2/ev_updated.txt")

conf_limit <- 1e-8

ev.fa <- ev.a %>%
  select(c('Sequence', 'PEP', 'pep_new', 'id', 'stan_peptide_id')) %>%
  filter(!is.na(pep_new)) %>%
  #filter(PEP > 0 & pep_new > 0 & PEP > conf_limit & pep_new > conf_limit) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

ev.fb <- ev.b %>%
  select(c('Sequence', 'PEP', 'pep_new', 'id', 'stan_peptide_id')) %>%
  filter(!is.na(pep_new)) %>%
  #filter(PEP > 0 & pep_new > 0 & PEP > conf_limit & pep_new > conf_limit) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

ev.fc <- ev.c %>%
  select(c('Sequence', 'PEP', 'pep_new', 'id', 'stan_peptide_id')) %>%
  filter(!is.na(pep_new)) %>%
  #filter(PEP > 0 & pep_new > 0 & PEP > conf_limit & pep_new > conf_limit) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

ev.fd <- ev.d %>%
  select(c('Sequence', 'PEP', 'pep_new', 'id', 'stan_peptide_id')) %>%
  filter(!is.na(pep_new)) %>%
  #filter(PEP > 0 & pep_new > 0 & PEP > conf_limit & pep_new > conf_limit) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

common_ids <- intersect(ev.fa$id, ev.fb$id)
common_ids <- intersect(common_ids, ev.fc$id)

ev.fa <- ev.fa %>% filter(id %in% common_ids) %>% arrange(id)
ev.fb <- ev.fb %>% filter(id %in% common_ids) %>% arrange(id)
ev.fc <- ev.fc %>% filter(id %in% common_ids) %>% arrange(id)
ev.fd <- ev.fd %>% filter(id %in% common_ids) %>% arrange(id)

# plot updated PEPs against each other ------------------------------------

k <- 100
contour_cols <- viridis::viridis(k)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

n <- sample.int(nrow(ev.fa), size=5e4)

#x <- log10(ev.fa$pep_new[n])
#y <- log10(ev.fb$pep_new[n])

x <- log10(ev.fa$pep_new[n])
y <- log10(ev.fc$pep_new[n])

dens <- get_density(x, y, k)

plot(x, y, pch=16, cex=0.25,
     # color by density
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
     
     xlim=c(-8, 0), ylim=c(-8, 0),
     #xlab='No Bootstrapping', ylab='Non-parametric Bootstrapping',
     xlab='Non-Parametric', ylab='Parametric',
     main='Updated PEP')
abline(a=0, b=1, col='red', lwd=2)


# plot delta PEPs from maxquant -------------------------------------------

set.seed(1)
n <- sample.int(nrow(ev.fa), 5e4)

x <- log10(ev.fa$PEP/ev.fa$pep_new)[n]
y <- log10(ev.fd$PEP/ev.fd$pep_new)[n]

inds <- ((x > 4) | (x < -4) | (y > 4) | (y < -4))
x <- x[-inds]
y <- y[-inds]

#k <- 200
#dens <- get_density(x, y, k)

#k <- 100
#cols <- colorRampPalette(c('blue', 'red'))

cols <- c(rgb(0,0,0,0.05), rgb(1,0,0,0.5))

plot(x, y, pch=16, cex=0.75,
     
     #col=rgb(0,0,0,0.05),
     #col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=100))],
     #col=paste0(cols(k)[findInterval(log10(ev.fa$PEP), seq(-5, 0, length.out=k))], '66'),
     col=cols[(ev.fa$PEP > 0.5)+1],
     
     xlim=c(-4, 4), ylim=c(-4, 4),
     #xlab='Parametric Bootstrap', 
     xlab='No Bootstrap', 
     ylab='Parametric Bootstrap\nMixture Sampling',
     main='Change in PEP\nlog10(PEP/Updated PEP)',
     mgp=c(2, 1, 0))
abline(a=0, b=1, col='red', lwd=2)
abline(h=0, lty=2)
abline(v=0, lty=2)


# PEP difference vs. MQ PEP -----------------------------------------------

n <- sample.int(nrow(ev.fa), 1e4)

x <- log10(ev.fa$PEP)[n]
y <- log10(ev.fb$pep_new / ev.fa$pep_new)[n]
col <- 

plot(x, y, pch=16, cex=0.25,
     col=rgb(0,0,0,0.1),
     xlim=c(-10, 0), ylim=c(-4, 4),
     xlab='MaxQuant PEP', ylab='Updated PEP\nBootstrap / No bootstrap',
     mgp=c(2, 0.75, 0))
abline(h=0, col='red', lwd=1)

text(-10, 2, 'Bootstrap Higher PEP', adj=c(0, 0.5))
text(-10, -2, 'Bootstrap Lower PEP', adj=c(0, 0.5))

# PEP difference vs. other factors ----------------------------------------

dpep <- log10(ev.fa$pep_new / ev.fb$pep_new)
# positive - parametric does better
# negative - parametric does worse
#hist(dpep)

num_obs <- ev.fa %>%
  filter(PEP < 0.5) %>%
  group_by(stan_peptide_id) %>%
  summarise(n=n()) %>%
  pull(n)

n_obs <- num_obs[ev.fa$stan_peptide_id+1]
n <- sample.int(nrow(ev.fa), 1e4)

plot(n_obs[n], dpep[n], pch=16, cex=0.5, col=rgb(0,0,0,0.1),
     xlim=c(0, 200), ylim=c(-4, 4),
     xlab='Number of samples', ylab='Change in Updated PEP',
     main='Non-parametric Bootstrap')
abline(a=0, b=0, col='red', lwd=2)
abline(v=0, col='black', lty=2)

text(50, 3, 'Upgrade')
text(50, -3, 'Downgrade')
