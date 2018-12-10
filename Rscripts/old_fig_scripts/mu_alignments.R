library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180730_1/ev_updated.txt')

peptide_params <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180730_1/peptide_params.txt')
exp_params <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180730_1/exp_params.txt')

ev.f <- ev %>%
  filter(!is.na(stan_peptide_id)) %>%
  filter(PEP < 0.5 & exclude=='False')

# apply linear regression params in reverse -------------------------------

ev.f <- ev.f %>%
  mutate(beta_0=exp_params$beta_0[exp_id+1],
         beta_1=exp_params$beta_1[exp_id+1],
         beta_2=exp_params$beta_2[exp_id+1],
         split_point=exp_params$split_point[exp_id+1])

#ev.f$jium <- ev.f$`Retention time`-ev.f$beta_0
ev.f$jium <- ev.f$`Retention time`
inds <- ev.f$mu < ev.f$split_point
ev.f$jium[inds] <- (ev.f$`Retention time`[inds]-ev.f$beta_0[inds]) / ev.f$beta_1[inds]

inds <- ev.f$mu > ev.f$split_point
ev.f$jium[inds] <- ((ev.f$`Retention time`[inds]-ev.f$beta_0[inds]-(ev.f$beta_1[inds]*ev.f$split_point[inds])) / ev.f$beta_2[inds]) + ev.f$split_point[inds]

plot(ev.f$mu[1:5e3], ev.f$jium[1:5e3])
abline(a=0, b=1, lwd=2, col='red')

# get mus from mean, median RTs -------------------------------------------

ev$mu

mus <- ev.f %>%
  group_by(stan_peptide_id) %>%
  mutate(weight=((1-PEP)-(1-0.5)) / 0.5) %>%
  summarise(mu_mean=mean(jium),
            mu_median=median(jium),
            mu_mean_weighted=sum(jium*weight) / sum(weight),
            mu_aligned=first(mu),
            num_exps=length(unique(`Raw file`)),
            num_obs=n()) %>%
  mutate(mu_mean_res=mu_mean, 
         mu_median_res=mu_median, 
         mu_mean_weighted_res=mu_mean_weighted,
         mu_init=peptide_params$init_mu) %>%
  mutate(mu_init_res=mu_init) %>%
  mutate_at(c('mu_mean_res', 'mu_median_res', 'mu_mean_weighted_res', 'mu_init_res'),
            funs(. - mu_aligned))



# residual histograms -----------------------------------------------------

layout(rbind(c(1, 2), c(3, 4)))

brks <- seq(-30, 30, by=0.05)

hist(mus$mu_mean_res, brks, freq=F,
     xlim=c(-2, 2), 
     xlab='Difference from Aligned Mu (min)',
     main='Mean Aligned RT')
hist(mus$mu_median_res, brks, freq=F,
     xlim=c(-2, 2), 
     xlab='Difference from Aligned Mu (min)',
     main='Median Aligned RT')
hist(mus$mu_mean_weighted_res, brks, freq=F,
     xlim=c(-2, 2), 
     xlab='Difference from Aligned Mu (min)',
     main='Mean Aligned RT (Weighted)')
hist(mus$mu_init_res, brks, freq=F,
     xlim=c(-2, 2), 
     xlab='Difference from Aligned Mu (min)',
     main='Mu -- Initial Value')

#dev.off()

# ridgeplot ---------------------------------------------------------------

mus$bin <- cut(mus$mu_aligned, seq(0, 60, by=5))

ggplot(mus) +
  geom_density_ridges(aes(x=mu_mean_res, y=bin, group=bin))
ggplot(mus) +
  geom_density_ridges(aes(x=mu_median_res, y=bin, group=bin))
ggplot(mus) +
  geom_density_ridges(aes(x=mu_init_res, y=bin, group=bin))

mus %>%
  gather(type, value, contains('res')) %>%
  ggplot() +
  geom_density_ridges(aes(x=value, y=type, group=type))
