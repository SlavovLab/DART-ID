library(scales)
library(readr)
library(tidyverse)
library("splines")
library(rstan)

evidence <- read_tsv("~/Google\ Drive/Ali_RT_bayesian/dat/ev.adj.txt")

## Filter of PEP  < .05 and remove REV and CONT and experiments on wrong LC column
subEvidence <- evidence %>% filter(PEP < 0.05) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>%
    filter(!grepl('REV*', `Leading razor protein`)) %>%
    filter(!grepl('CON*',`Leading razor protein`))  %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP")

## Add factor indices
subEvidence <- subEvidence %>% mutate(exp_id=`Raw file`) %>% mutate_at("exp_id", funs(as.numeric(as.factor(.))))

## Look at mean retention times to get a sense of alignment
rt_means <- subEvidence %>% group_by(`Peptide ID`, `exp_id`) %>% summarise(avg=median(`Retention time`)) ##, pep_avg=mean(`PEP`))
rt_means <- rt_means %>% spread(key=`exp_id`, value=avg)


rank_means <- subEvidence %>% group_by(`exp_id`) %>% mutate_at(vars(`Retention time`), funs(rank(., na.last="keep")/(n()+1))) %>% ungroup() %>% group_by(`exp_id`, `Peptide ID`) %>% summarise(avg=median(`Retention time`)) %>% spread(key=`exp_id`, value=avg)

num_experiments <- length(unique(subEvidence[["exp_id"]]))

rt_means %>%  arrange(`33`) %>% ggplot(aes(x=`33`, y=`43`)) + geom_point() + geom_abline(a=0, b=1)
rt_means %>%  ggplot(aes(x=`38`, y=`33`)) + geom_line() + geom_abline(a=0, b=1)
rt_means  %>% ggplot(aes(x=`99`, y=`68`)) + geom_point() + geom_abline(a=0, b=1)
rt_means  %>% ggplot(aes(x=`13`, y=`44`)) + geom_point() + geom_abline(a=0, b=1)


rank_means %>%  ggplot(aes(x=`38`, y=`33`)) + geom_point() + geom_abline(a=0, b=1)
rank_means  %>% ggplot(aes(x=`99`, y=`68`)) + geom_point() + geom_abline(a=0, b=1)
rank_means  %>% ggplot(aes(x=`13`, y=`44`)) + geom_point() + geom_abline(a=0, b=1)
rank_means  %>% ggplot(aes(x=`50`, y=`60`)) + geom_point() + geom_abline(a=0, b=1)

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
muij_map <- match(paste(stan_peptide_id, exp_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

data <- list(num_experiments=num_experiments, num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, muij_map=muij_map,
             muij_to_exp=muij_to_exp, muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=exp_id, peptide_id=stan_peptide_id,
             retention_times=retention_times)

sm <- stan_model(file="fit_RT2.stan")

muInit <- sapply(unique(stan_peptide_id), function(pid) mean(retention_times[stan_peptide_id==pid]))
initList <- list(mu=muInit,
                 beta_0 = rep(1, num_experiments),
                 beta_1 = rep(1, num_experiments),
                 beta_2 = rep(1, num_experiments),
                 sigma_slope=rep(0.1, num_experiments),
                 sigma_intercept=rep(3, num_experiments),
                 split_point=rep(median(muInit), num_experiments))

                 

pars <- optimizing(sm, data=data, init=initList, iter=10000, verbose=TRUE)$par

## pars <- sampling(sm, data=data, init=list(initList), chains=1, verbose=TRUE)$par
save(pars, file="params_rt2.RData")

sm2 <- stan_model(file="fit_RT_mu_fixed.stan")
data2 <- data
data2$mu <- pars[grep("mu\\[", names(pars))]
pars <- optimizing(sm2, data=data2, init=initList, iter=10000, verbose=TRUE)$par

## Beta fit
muij_fit <- pars[grep("muij", names(pars))]
for(exp in 1:num_experiments) {
    c(pars[paste0("beta_1[", exp, "]")], pars[paste0("beta_2[", exp, "]")])

    exp_indices <-  muij_to_exp==exp
    pep_indices <- muij_to_pep[exp_indices]
    split <-   pars[grep(paste0("split_point\\[", exp, "\\]"), names(pars))]

    exp_fit <- muij_fit[exp_indices]
    sigmas <- pars["sigma_intercept[20]"] + pars["sigma_slope[20]"]/100 * sort(exp_fit)
    min(sigmas)
    max(sigmas)


    predicted <- (muij_fit[muij_map])[(muij_to_exp==exp)[muij_map]]
    observed <- retention_times[(muij_to_exp==exp)[muij_map]]
    residual <- observed - predicted

    pdf(sprintf("figs/second-fit-%i.pdf", exp))
    plot(predicted, observed, pch=19, cex=0.2)
    abline(a=0, b=1, col="blue")
    abline(v=split, col="blue")
    dev.off()

}
