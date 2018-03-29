library(scales)
library(readr)
library(tidyverse)
library(splines)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

evidence <- read_tsv('dat/evidence_elite.txt')

evidence <- evidence %>%
  rename(`Sequence ID`=`Peptide ID`) %>%
  rename(`Peptide ID`=`Mod. peptide ID`) # alias the modified peptide ID as the peptide ID

## Filter of PEP < pep_thresh and remove REV and CON and experiments on wrong LC column
pep_thresh <- 0.5

# load exclusion list - common contaminants + keratins
exclude <- read_lines('pd_exclude.txt')

subEvidence <- evidence %>% 
  # remove EB/ES experiments - human only
  filter(!grepl('28C|28D|30[K-L]|31[A-F]|32[E-I]|3[3-6][A-I]', `Raw file`)) %>%
  filter(PEP < pep_thresh) %>%
  filter(!grepl('REV*', `Leading razor protein`)) %>%
  filter(!grepl('CON*',`Leading razor protein`))  %>%
  filter(!grepl(paste(exclude, collapse='|'), Proteins)) %>%
  filter(`Retention length` < 5) %>%
  select("Peptide ID", "Raw file", "Retention time", "PEP") %>%
  ## Add factor indices
  mutate(exp_id=as.numeric(as.factor(`Raw file`)))

## Remove peptides that occur in less than n experiments
toRemove <- subEvidence %>% 
  group_by(`Peptide ID`) %>% 
  summarise(num.exps=length(unique(`exp_id`))) %>% 
  filter(num.exps <= 10) %>%
  pull(`Peptide ID`)
subEvidence <- subEvidence %>% 
  filter(!(`Peptide ID` %in% toRemove))

num_experiments <- max(subEvidence[["exp_id"]])

## "true peptide id" matches peptide id in evidence file
## "peptide id" is index in 1:num_peptides for stan
raw_peptide_id <- subEvidence[["Peptide ID"]]
stan_peptide_id <- as.numeric(as.factor(subEvidence[["Peptide ID"]]))
subEvidence[["Stan ID"]] <- stan_peptide_id

num_total_observations <- nrow(subEvidence)
num_peptides <- length(unique(stan_peptide_id))

retention_times <- subEvidence[["Retention time"]]
exp_id <- subEvidence[["exp_id"]]

pep_exp_all <- paste(stan_peptide_id, exp_id, sep=" - ")
pep_exp_pairs <- unique(pep_exp_all)
num_pep_exp_pairs <- length(pep_exp_pairs)

## map peptides to stan peptides in stan
muij_map <- match(paste(stan_peptide_id, exp_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

## Remove zeros
pep <- subEvidence$PEP
pep <- pep + .Machine$double.eps
#pep <- pep + 1e-8

## Gather all data required by the stan model into a list
data <- list(num_experiments=num_experiments, 
             num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, 
             muij_map=muij_map,
             muij_to_exp=muij_to_exp, 
             muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=exp_id, 
             peptide_id=stan_peptide_id,
             retention_times=retention_times,
             mean_log_rt=mean(log(retention_times)),
             sd_log_rt=sd(log(retention_times)),
             pep = pep,
             max_retention_time=max(retention_times))

# fix randomness
set.seed(0)

# Initialize canonical retention time priors for each peptide
muInit <- sapply(unique(stan_peptide_id), function(pid) {
  # mean of peptide retention times, weighted by PEP
  weights <- ((1 - pep[stan_peptide_id == pid]) - (1 - pep_thresh)) / pep_thresh
  sum(retention_times[stan_peptide_id==pid] * weights) / sum(weights) + rnorm(1, 0, 10)
})
# negative or very low retention times not allowed. floor at 5 minutes
mu.min <- 5
muInit[muInit <= mu.min] <- mu.min
# canonical retention time shouldn't be bigger than largest real RT
mu.max <- max(retention_times)
muInit[muInit > mu.max] <- mu.max

# take retention times and distort by +- 10 mins
rt_distorted <- retention_times + rnorm(num_total_observations, 0, 10)
# make sure distorted retention times stay within bounds of real ones
rt_distorted[rt_distorted > max(retention_times)] <- max(retention_times)
rt_distorted[rt_distorted < min(retention_times)] <- min(retention_times)

# initialize priors for the linear regression
# first element of vector is beta_0, or the intercept
# second element is beta_1, the slope
beta_init <- rbind(rep(10, num_experiments), rep(1, num_experiments))

for( i in 1:10 ) {
  print(i) 
  
  # for each experiment, fit a simple linear regression
  # between the distorted RTs and the initial canonical retention times
  lm_coefs <- sapply(sort(unique(exp_id)), function(id) {
    rt_cur <- rt_distorted[exp_id == id] 
    mu_cur <- muInit[match(stan_peptide_id[exp_id == id], unique(stan_peptide_id))]
    pep_cur <- pep[exp_id == id]
    lm_cur <- lm(rt_cur ~ mu_cur, weights = 1-pep_cur)
    
    lm_cur$coefficients
  })
  # store linear regression parameters
  beta_init[1, ] <- lm_coefs[1, ]
  beta_init[2, ]  = lm_coefs[2, ]
  # store this set of priors before iterating again
  muInitPrev <- muInit
  
  # calculate new set of canonical RTs based on linear regression params
  mu_pred <- (rt_distorted - beta_init[1, exp_id]) / beta_init[2, exp_id] 
  # make sure new canonical RTs are within same range as distorted RTs
  mu_pred[mu_pred <= 0] <- min(rt_distorted)
  mu_pred[mu_pred >= max(rt_distorted)] <- max(rt_distorted)
  
  muPrev <- muInit
  
  # new set of priors for canonical RTs based on weighted combination of
  # this set of predicted canonical RTs
  muInit <- sapply(unique(stan_peptide_id), function(pid) {
    weights <- ((1-pep[stan_peptide_id==pid]) - (1-pep_thresh))/pep_thresh
    sum(mu_pred[stan_peptide_id==pid] * weights) / sum(weights)
  })
  
  print(sum(muPrev - muInit)^2 / length(muInit))
}

beta_0 <- beta_init[1, ]
beta_1 <- beta_init[2, ]
# apply lower bound of (-1.0 * min(beta_1) * min(muInit)) to beta_0
# where (-1.0 * min(beta_1) * min(muInit)) is the lowest possible intercept
# given the lowest possible mu and lowest possible beta_1
beta_0[beta_0 <= -1.0 * min(beta_1) * min(muInit)] <- (-1.0 * min(beta_1) * min(muInit)) + 1e-3

# apply upper bound to prior canonical RTs
muInit[muInit >= max(retention_times)] <- 0.95 * max(retention_times)

# create prior list for STAN
initList <- list(mu=muInit,
                 beta_0 = beta_0,
                 beta_1 = beta_1,
                 sigma_slope=rep(0.1, num_experiments),
                 sigma_intercept=rep(0.1, num_experiments))

save(data, initList, file='dat/alignment_data_linear.RData')

## Compile the stan code
#sm <- stan_model(file="fit_RT3.stan")
sm <- stan_model(file='fit_RT_linear.stan')
save(sm, file='fit_RT_linear.RData')

#load('dat/fit_RT3.RData')
#load('fit_RT_linear.RData')
#load('dat/alignment_data_linear.RData')

start <- Sys.time()
pars <- optimizing(sm, data=data, init=initList, iter=20000, verbose=TRUE)$par
print(Sys.time() - start)

save(pars, file='dat/params.Fit.linear.RData')
#save(pars, file='dat/params.Fit3c.RData')

## Beta fit
library(rmutil)
muij_fit <- pars[grep("muij", names(pars))]
sigma_ij <- pars[grep("sigma_ij", names(pars))]

residual_vec <- c()
exp_vec <- c()
col_vec <- c()

pep_col_code <- cut(pep, breaks=10)

for(exp in 1:num_experiments) {
  betas <- c(pars[paste0("beta_0[", exp, "]")], 
             pars[paste0("beta_1[", exp, "]")], 
             pars[paste0("beta_2[", exp, "]")])
  
  ## Get canonical retention times
  mu <-  pars[grep("mu\\[", names(pars))]
  
  ## get experiment indices
  exp_indices <- muij_to_exp==exp
  
  ## split points
  split <-  pars[grep(paste0("split_point\\[", exp, "\\]"), names(pars))]
  
  ## estimated muij for experiment
  exp_fit <- muij_fit[exp_indices]
  sigmas <- pars["sigma_intercept[20]"] + pars["sigma_slope[20]"]/100 * sort(exp_fit)
  #min(sigmas)
  #max(sigmas)
  
  ## predicted values (duplicated because multiple observations)
  predicted <- (muij_fit[muij_map])[exp_indices[muij_map]]
  ## predicted sd (duplicated because multiple observations)
  predicted_sd <- (sigma_ij[muij_map])[(muij_to_exp==exp)[muij_map]]
  
  mus <- mu[(muij_to_pep[muij_map])[exp_indices[muij_map]]]
  observed <- retention_times[exp_indices[muij_map]]
  obs_peps <- pep[exp_indices[muij_map]]
  obs_code <- pep_col_code[exp_indices[muij_map]]
  residual <- observed - predicted
  
  residual_vec <- c(residual_vec, residual)
  col_vec <- c(col_vec, obs_code)
  exp_vec <- c(exp_vec, rep(exp, length(residual)))
  
  pdf(sprintf("tmp_figs/second-fit-%i.pdf", exp))
  plot(predicted, observed, pch=19, cex=0.2)
  abline(a=0, b=1, col="blue")
  abline(v=split, col="blue")
  dev.off()
  
  
  pdf(sprintf("tmp_figs/canonical_v_original-%i.pdf", exp))
  plot(mus, observed, pch=19, cex=0.2)
  abline(v=split, col="blue", lty=2)
  segments(x0=0, y0=betas[1], x1=split, y1=betas[1]+betas[2]*split, col="green", lwd=1.5)
  segments(x0=split, y0=betas[1]+betas[2]*split, x1=400, y1=betas[1] + betas[2]*split +betas[3]*(400-split), col="red", lwd=1.5)
  dev.off()
  
  cols <- rev(heat.colors(10))
  
  pdf(sprintf("tmp_figs/residuals-fit-%i.pdf", exp))
  plot(predicted, observed-predicted,pch=19, cex=0.2, col=cols[obs_code])
  lines(predicted[order(predicted)],
        sapply(predicted_sd, function(s) qlaplace(.025, 0, s))[order(predicted)],
        col="red")
  
  lines(predicted[order(predicted)],
        sapply(predicted_sd, function(s) qlaplace(.975, 0, s))[order(predicted)],
        col="red")
  
  
  abline(v=split, col="blue")
  dev.off()
  
}

library(ggridges)
df <- data.frame(x=abs(residual_vec), y=as.factor(col_vec))

df %>% group_by(y) %>% summarise(mean=mean(x), qtl=quantile(x, 0.99))
head(df)
pdf("~/Desktop/tmp.pdf")
ggplot(df, aes(x=log2(x), y=y, fill=y)) + geom_density_ridges()
dev.off() 
summary(pars[grep("sigma_intercept", names(pars))])
summary(rlnorm(1000, 0, 2))


subEvidence %>% filter(exp_id == 50) %>% filter(duplicated(`Peptide ID`, fromLast=TRUE))
tail(pep_col_code)
pep[pep_col_code==5]

summary(retention_times)
