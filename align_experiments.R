##################################################################
## Fit piecewise linear model for aligning retention times
##################################################################

library(rstan)
library(scales)
library(readr)

evidence <- read_tsv("~/Google Drive/Ali_RT_bayesian/dat/evidence.txt")
subEvidence <- evidence[evidence$PEP < 0.05, c("Peptide ID", "Raw file", "Retention time")]

experimentFactors <- as.factor(subEvidence[["Raw file"]])
experiment_id <- as.numeric(experimentFactors)
num_experiments <- max(experiment_id)

## WARNING:
## "true_peptide_id" matches peptide id in evidence file
## "peptide_id" is index in 1:num_peptides which is the input for stan
## ultimately need to convert between the two but below we work with peptide_id

true_peptide_id <- subEvidence[["Peptide ID"]]
peptide_id <- as.numeric(as.factor(subEvidence[["Peptide ID"]]))

num_total_observations <- nrow(subEvidence)
num_peptides <- length(unique(peptide_id))

retention_times <- subEvidence[["Retention time"]]

pep_exp_all <- paste(peptide_id, experiment_id, sep=" - ")
pep_exp_pairs <- unique(pep_exp_all)
num_pep_exp_pairs <- length(pep_exp_pairs)
muij_map <- match(paste(peptide_id, experiment_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

## Gather all data required by the stan model into a list
data <- list(num_experiments=num_experiments, num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, muij_map=muij_map,
             muij_to_exp=muij_to_exp, muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=experiment_id, peptide_id=peptide_id,
             retention_times=retention_times)

## Compile the stan code
sm <- stan_model(file="fit_RT.stan")

## Initialize the parameters to something sensible
muInit <- aggregate(`Retention time` ~ `Peptide ID`, 
                    data=subEvidence, FUN=mean)$`Retention time`
splitInit <- rep(mean(muInit), num_experiments)
initList <- list(mu=muInit,
                 beta_0=rep(0, num_experiments),
                 beta_1=rep(1, num_experiments),
                 beta_2=rep(1, num_experiments),
                 split_point=splitInit,
                 sigma_global=1, sigma=rep(1, num_peptides))

## Optimze and save params
pars <- optimizing(sm, data=data, init=initList, iter=20000, verbose=TRUE)$par
save(pars, file="params.RData")

## take a protein id and return the predicted mean retention time for all experiments
getPredicted <- function(pid) {
    exps <- experiment_id[peptide_id == pid]

    beta_0 <- pars[sprintf("beta_0[%i]", exps)]
    beta_1 <- pars[sprintf("beta_1[%i]", exps)]
    beta_2 <- pars[sprintf("beta_2[%i]", exps)]
    split_point <- pars[sprintf("split_point[%i]", exps)]
        
    mu <- pars[sprintf("mu[%s]", pid)]

    beta_0 + ifelse(mu < split_point,
                    beta_1*mu,
                    beta_1*split_point + beta_2*(mu - split_point))
}

## plot alignment of up to 4 experiments
checkExperiment <- function(experiments) {
    experiments <- experiments[1:min(length(experiments), 4)]
    quartz()
    par(mfrow=c(2, 4))
    count <- 0
    for(exp in experiments) {

        par(mfg=c(count %% 2 + 1, floor(count / 2) + 1))

        beta_0 <- pars[sprintf("beta_0[%i]", exp)]
        beta_1 <- pars[sprintf("beta_1[%i]", exp)]
        beta_2 <- pars[sprintf("beta_2[%i]", exp)]
        split_point <- pars[sprintf("split_point[%i]", exp)]
        
        meanRT <- sapply(unique(peptide_id[experiment_id == exp]), function(id) {
            mean(retention_times[experiment_id == exp & peptide_id == id])
        })
        names(meanRT) <- unique(peptide_id[experiment_id == exp])
        
        mus <- pars[sprintf("mu[%i]", unique(peptide_id[experiment_id == exp]))]
        plot(meanRT, beta_0[1] + ifelse(mus < split_point,
                beta_1*mus, beta_1*split_point + beta_2*(mus - split_point)), 
                main=sprintf("2 piece linear (%i)", exp), 
                ylab="predicted")
        abline(a=0, b=1, col="red")
        count <- count + 1
        
        par(mfg=c(count %% 2 + 1, floor(count / 2) + 1))
        plot(meanRT, mus, main=sprintf("Linear (%i)", exp), ylab="raw")
        coef <- lm(mus ~ meanRT)$coefficients
        abline(a=coef[1], b=coef[2], col="red")
        count <- count + 1
    }
}

checkExperiment(sample(140, 4))

checkExperiment(c(80, 6, 62, 55))

## Check means and variances for some experiments
peptideCounts <- sort(table(peptide_id), decreasing=TRUE)
par(mfrow=c(2, 4))
for(pid in names(peptideCounts)[101:108]) {
    rts <- retention_times[peptide_id == pid]
    diff <- rts - getPredicted(pid)
    hist(diff, breaks=50, freq=FALSE, main=pid)
    sigma <- pars[sprintf("sigma[%s]", pid)]
    curve(dnorm(x, 0, sigma),
          from = min(diff) - 0.1 * range(diff),
          max(diff) + 0.1 * range(diff), add=TRUE, col="red")
     print(sd(diff))

}




