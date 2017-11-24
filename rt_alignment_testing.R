library(scales)
library(readr)
library(tidyverse)
library(splines)
library(rstan)

## Old
## evidence <- read_tsv("~/Google\ Drive/Ali_RT_bayesian/dat/ev.adj.txt")

evidence <- read_tsv("~/Google\ Drive/bayesian_RT/dat/evidence_elite.txt")

## Filter of PEP  < .05 and remove REV and CONT and experiments on wrong LC column

pep_thresh <- 0.5

subEvidence <- evidence %>% filter(PEP < pep_thresh) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>%
    filter(!grepl('REV*', `Leading razor protein`)) %>%
    filter(!grepl('CON*',`Leading razor protein`))  %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP")



## Add factor indices
subEvidence <- subEvidence %>% mutate(exp_id=`Raw file`) %>% mutate_at("exp_id", funs(as.numeric(as.factor(.))))

## Remove peptides that occur in only one experiment
toRemove <- subEvidence %>% group_by(`Peptide ID`) %>% summarise(n=length(unique(`exp_id`))) %>% filter(n <= 10)
subEvidence <- subEvidence %>% filter(!(`Peptide ID` %in% toRemove[["Peptide ID"]]))

## Look at mean retention times to get a sense of alignment
rt_means <- subEvidence %>% group_by(`Peptide ID`, `exp_id`) %>% summarise(avg=median(`Retention time`)) ##, pep_avg=mean(`PEP`))
rt_means <- rt_means %>% spread(key=`exp_id`, value=avg)


rank_means <- subEvidence %>% group_by(`exp_id`) %>% mutate_at(vars(`Retention time`), funs(rank(., na.last="keep")/(n()+1))) %>% ungroup() %>% group_by(`exp_id`, `Peptide ID`) %>% summarise(avg=median(`Retention time`)) %>% spread(key=`exp_id`, value=avg)

num_experiments <- length(unique(subEvidence[["exp_id"]]))

rt_means %>%  arrange(`33`) %>% ggplot(aes(x=`33`, y=`43`)) + geom_point() + geom_abline()
rt_means %>%  ggplot(aes(x=`38`, y=`33`)) + geom_line() + geom_abline()
rt_means  %>% ggplot(aes(x=`99`, y=`68`)) + geom_point() + geom_abline()
rt_means  %>% ggplot(aes(x=`13`, y=`44`)) + geom_point() + geom_abline()
rt_means  %>% ggplot(aes(x=`50`, y=`60`)) + geom_point() + geom_abline()


rank_means %>%  ggplot(aes(x=`38`, y=`33`)) + geom_point() + geom_abline()
rank_means  %>% ggplot(aes(x=`99`, y=`68`)) + geom_point() + geom_abline()
rank_means  %>% ggplot(aes(x=`13`, y=`44`)) + geom_point() + geom_abline()
rank_means  %>% ggplot(aes(x=`50`, y=`60`)) + geom_point() + geom_abline()

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
pep <- pep + 1e-8

## Gather all data required by the stan model into a list
data <- list(num_experiments=num_experiments, num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, muij_map=muij_map,
             muij_to_exp=muij_to_exp, muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=exp_id, peptide_id=stan_peptide_id,
             retention_times=retention_times,
             pep = pep,
             max_retention_time=max(retention_times))

## Compile the stan code
sm <- stan_model(file="fit_RT3.stan")

muInit <- sapply(unique(stan_peptide_id), function(pid) {
    weights <- ((1-pep[stan_peptide_id==pid]) - (1-pep_thresh))/pep_thresh
    sum(retention_times[stan_peptide_id==pid] * weights) / sum(weights) + rnorm(1, 0, 10)
})

muInit[muInit <= 0] <- 5

max_mu <- max(retention_times)
muInit[muInit > max_mu] <- max_mu

rt_distorted <- retention_times + rnorm(num_total_observations, 0, 10)
rt_distorted[rt_distorted > max(retention_times)] <- max(retention_times)
rt_distorted[rt_distorted < min(retention_times)] <- min(retention_times)
beta_init <-  rbind(rep(10, num_experiments), rep(1, num_experiments))

for( i in 1:10 ) {

    print(i) 

    lm_coefs <- sapply(sort(unique(exp_id)), function(id) {
        rt_cur <- rt_distorted[exp_id == id] 
        mu_cur <- muInit[match(stan_peptide_id[exp_id == id], unique(stan_peptide_id))]
        pep_cur <- pep[exp_id == id]
        lm_cur <- lm(rt_cur ~ mu_cur, weights = 1-pep_cur)

        lm_cur$coefficients
    })
    
    beta_init[1, ] <- lm_coefs[1, ]
    beta_init[2, ]  = lm_coefs[2, ]
    
    muInitPrev <- muInit

    mu_pred <- (rt_distorted - beta_init[1, exp_id]) / beta_init[2, exp_id] 
    
    mu_pred[mu_pred <= 0] <- min(rt_distorted)
    mu_pred[mu_pred >= max(rt_distorted)] <- max(rt_distorted)
    muPrev <- muInit
    muInit <- sapply(unique(stan_peptide_id), function(pid) {
        weights <- ((1-pep[stan_peptide_id==pid]) - (1-pep_thresh))/pep_thresh
        sum(mu_pred[stan_peptide_id==pid] * weights) / sum(weights)
    })

    print( sum(muPrev - muInit)^2/length(muInit))
##    plot(beta_init[1, emuInit[stan_peptide_id], rt_distorted, pch=19, cex=0.1)
##    id <- 20; plot(muInit[stan_peptide_id[exp_id==id]]*beta_init[2, id] + beta_init[1, id], rt_distorted[exp_id == id], pch=19, cex=0.1)

}

beta_0 = beta_init[1, ]
beta_1 <- beta_2 <- beta_init[2, ]
beta_0[beta_0 <= -1.0*min(beta_1)*min(muInit)] <- -1.0*min(beta_1)*min(muInit) + 1e-3
       
muInit[muInit >= max(retention_times)] <- 0.95*max(retention_times)


initList <- list(mu=muInit,
                 beta_0 = beta_0,
                 beta_1 = beta_1,
                 beta_2 = beta_2,
                 sigma_slope=rep(0.1, num_experiments),
                 sigma_intercept=rep(0.1, num_experiments),
                 split_point=rep(median(muInit), num_experiments))

start <- Sys.time()
pars <- optimizing(sm, data=data, init=initList, iter=20000, verbose=TRUE)$par
print(Sys.time() - start)

## Beta fit
library(rmutil)
muij_fit <- pars[grep("muij", names(pars))]
sigma_ij <- pars[grep("sigma_ij", names(pars))]

residual_vec <- c()
exp_vec <- c()
col_vec <- c()

pep_col_code <- cut(pep, breaks=10)

for(exp in 1:num_experiments) {
    betas <- c(pars[paste0("beta_0[", exp, "]")], pars[paste0("beta_1[", exp, "]")], pars[paste0("beta_2[", exp, "]")])

    ## Get canonical retention times
    mu <-  pars[grep("mu\\[", names(pars))]

    ## get experiment indices
    exp_indices <-  muij_to_exp==exp

    ## split points
    split <-  pars[grep(paste0("split_point\\[", exp, "\\]"), names(pars))]

    ## estimated muij for experiment
    exp_fit <- muij_fit[exp_indices]
    sigmas <- pars["sigma_intercept[20]"] + pars["sigma_slope[20]"]/100 * sort(exp_fit)
    min(sigmas)
    max(sigmas)

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
