library(scales)
library(readr)
library(tidyverse)
library("splines")
library(rstan)

evidence <- read_tsv("~/Google\ Drive/Ali_RT_bayesian/dat/ev.adj.txt")

hist(evidence[["Retention time"]], breaks=50, probability=TRUE)
abline(v=mean(evidence[evidence[["Peptide ID"]] == 208, ][["Retention time"]]),
       lwd=3, col="blue")
abline(v=mean(evidence[evidence[["Peptide ID"]] == 290, ][["Retention time"]]),
       lwd=3, col="red")
abline(v=mean(evidence[evidence[["Peptide ID"]] == 196, ][["Retention time"]]),
       lwd=3, col="green")


which.max(table(evidence[["Leading razor protein"]]))

sub <- evidence[evidence[["Leading razor protein"]] == "sp|P36578|RL4_HUMAN", ]

######################
## Peptide Exploration
######################

peptideCounts <- sort(table(evidence["Peptide ID"]), decreasing=TRUE)
hist(log10(as.numeric(peptideCounts)))
head(peptideCounts, n=100)

## to look at
## 592348, 572738

## explore some peptides
pep <- 592348
fdr <- pmin(evidence$PEP[evidence["Peptide ID"] == pep], 1)
print(fdr)

## hist(evidence[["Retention time"]], breaks=50, probability=TRUE, ylim=c(0, 0.01))
rts <- evidence[evidence[["Peptide ID"]] == pep, ][["Retention time"]]
hist(evidence[["Retention time"]], breaks=50, probability=TRUE)
abline(v=rts, col=alpha("blue", alpha=(1-fdr)))

######################
## Protein Exploration
######################

proteinCounts <- sort(table(evidence["Leading razor protein"]), decreasing=TRUE)
## for(protein in names(proteinCounts)) {
protein <- "sp|P19338|NUCL_HUMAN"
proteinEvidence <- evidence[evidence["Leading razor protein"] == protein, ]
table(proteinEvidence[["Peptide ID"]])

## 474, 480 , 467437
i1 <- proteinEvidence[proteinEvidence["Peptide ID"] == 474, c("Raw file", "Intensity", "PEP")]
i2 <- proteinEvidence[proteinEvidence["Peptide ID"] == 418167, c("Raw file", "Intensity", "PEP")]

indices <- intersect(i1[["Raw file"]], i2[["Raw file"]])
res1 <- log(i1[match(indices, i1[["Raw file"]]), "Intensity"])
res2 <- log(i2[match(indices, i2[["Raw file"]]), "Intensity"])

pep1 <- pmin(i1[match(indices, i1[["Raw file"]]), ]$PEP, 1)
pep2 <- pmin(i2[match(indices, i2[["Raw file"]]), ]$PEP, 1)
sort(pep1)
sort(pep2)
cor.test(pep1, pep2)

cor.test(res1$Intensity[-1], res2$Intensity[-1], pch=19)
plot(res1$Intensity, res2$Intensity, pch=19)
plot(res1$Intensity, res2$Intensity, pch=19, col=alpha("black", alpha=abs(pep1-pep2)))


##}

## for( rf in unique(evidence[["Raw file"]])) {
##     print(summary(log(evidence$Intensity[evidence[["Raw file"]] == rf]), na.rm=TRUE))
## }


######################
## Test STAN
######################

## Filter of PEP  < .05 and remove REV and CONT and experiments on wrong LC column
subEvidence <- evidence %>% filter(PEP < 0.05) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>%
    filter(!grepl('REV*', `Leading razor protein`)) %>%
    filter(!grepl('CON*',`Leading razor protein`))  %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP")

## Add factor indices
subEvidence <- subEvidence %>% mutate(exp_id=`Raw file`) %>% mutate_at("exp_id", funs(as.numeric(as.factor(.))))


## Look at mean retention times to get a sense of alignment
rt_means <- subEvidence %>% group_by(`Peptide ID`, `exp_id`) %>% summarise(avg=mean(`Retention time`)) ##, pep_avg=mean(`PEP`))
rt_means <- rt_means %>% spread(key=`exp_id`, value=avg)
num_experiments <- length(unique(subEvidence[["exp_id"]]))

rt_means %>%  arrange(`33`) %>% ggplot(aes(x=`33`, y=`43`)) + geom_line()
rt_means %>%  ggplot(aes(x=`38`, y=`33`)) + geom_line()
rt_means  %>% ggplot(aes(x=`99`, y=`68`)) + geom_line()
rt_means  %>% ggplot(aes(x=`13`, y=`44`)) + geom_line()

num_experiments <- max(subEvidence[["exp_id"]])

## "true peptide id" matches peptide id in evidence file
## "peptide id" is index in 1:num_peptides for stan
true_peptide_id <- subEvidence[["Peptide ID"]]
peptide_id <- as.numeric(as.factor(subEvidence[["Peptide ID"]]))

num_total_observations <- nrow(subEvidence)
num_peptides <- length(unique(peptide_id))

retention_times <- subEvidence[["Retention time"]]
exp_id <- subEvidence[["exp_id"]]

pep_exp_all <- paste(peptide_id, exp_id, sep=" - ")
pep_exp_pairs <- unique(pep_exp_all)
num_pep_exp_pairs <- length(pep_exp_pairs)
muij_map <- match(paste(peptide_id, exp_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

## For basis splines
B <- t(bs(X, knots=seq(0, max(retention_times), length.out=5),
          degree=3, intercept = TRUE))
sm<-stan_model("fit_basis.stan")
fit<-sampling(sm,iter=500,control = list(adapt_delta=0.95))

data <- list(num_experiments=num_experiments, num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, muij_map=muij_map,
             muij_to_exp=muij_to_exp, muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=exp_id, peptide_id=peptide_id,
             retention_times=retention_times)

## sm <- stan_model(file="fit_RT.stan")

## sm <- stan_model(file="fit_RT_weighted.stan")
## sm <- stan_model(file="fit_RT_laplace.stan")
sm <- stan_model(file="fit_RT_laplace_betacdf.stan")
## sm <- stan_model(file="fit_RT_laplace_splines.stan")


data <- list(num_experiments=num_experiments, num_peptides=num_peptides,
             num_pep_exp_pairs=num_pep_exp_pairs, muij_map=muij_map,
             muij_to_exp=muij_to_exp, muij_to_pep=muij_to_pep,
             num_total_observations=num_total_observations,
             experiment_id=exp_id, peptide_id=peptide_id,
             retention_times=retention_times)



muInit <- sapply(unique(peptide_id), function(pid) mean(retention_times[peptide_id==pid]))
muInit <- rep(mean(retention_times), length(unique(peptide_id)))
mu_beta <-  rep(1/2, num_experiments)
aplusb <-  rep(2, num_experiments)
gammaInit <- rep(200, num_experiments)
initList <- list(mu=muInit/250,
                 log_mu_beta=log(mu_beta),
                 log_aplusb=log(aplusb),
                 gamma=gammaInit,
                 intercept=rep(0, num_experiments), 
                 sigma_slope=rep(1, num_experiments),
                 sigma_slope_global=1,
                 sigma_ij=rep(1, num_pep_exp_pairs))
                 

pars <- optimizing(sm, data=data, init=initList, iter=1000, verbose=TRUE)$par

pars <- sampling(sm, data=data, init=list(initList), chains=1, verbose=TRUE)$par
save(pars, file="params_laplace.RData")





## Beta fit
muij_fit <- pars[grep("muij", names(pars))]
exp <- 2
exp_indices <- muij_to_pep[muij_to_exp==exp]
exp_fit <- muij_fit[exp_indices]

muInit[muij_to_pep]

meanRT <- subEvidence %>% filter(exp_id==exp) %>% 
    group_by(`Peptide ID`) %>%
    summarise(avg=mean(`Retention time`)) %>%
    ungroup() %>% 
    select("avg") %>% unlist

plot(exp_fit, meanRT)



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

checkExperiment <- function(experiments, pdf_file=NULL, residual=TRUE, acf=FALSE, width=14) {

    if(!is.null(pdf)) {
        pdf(pdf_file, width=width)
        print(pdf_file)
    }

    par(mfrow=c(2, 4))
    count <- 0

    for(exp in experiments) {

        print("---------")
        print(floor(count / 2) + 1)
        par(mfg=c(count %% 2 + 1, floor(count / 2) + 1))

        beta_0 <- pars[sprintf("beta_0[%i]", exp)]
        beta_1 <- pars[sprintf("beta_1[%i]", exp)]
        beta_2 <- pars[sprintf("beta_2[%i]", exp)]
        split_point <- pars[sprintf("split_point[%i]", exp)]

        experimentEvidence <- subEvidence %>% filter(`exp_id` == exp)
        
        meanRT <- experimentEvidence %>%
            group_by(`Peptide ID`, `exp_id`) %>%
            summarise(avg=mean(`Retention time`)) %>%
            ungroup() %>% 
            select("avg") %>% unlist

        experiment_id <- subEvidence[["exp_id"]]
        mus <- pars[sprintf("mu[%i]", unique(peptide_id[experiment_id == exp]))]
        predicted <- beta_0[1] + ifelse(mus < split_point, 
                                        beta_1*mus, beta_1*split_point + beta_2*(mus - split_point))
        browser()
        if(residual) {
            ord <- order(meanRT)
            plot(meanRT[ord], predicted[ord] - meanRT[ord],  main=sprintf("2 piece linear (%i)", exp), ylab="predicted", type="l")
            abline(h=0, b=1, col="red", lty=2)
        } else if(acf) {
            ord <- order(meanRT)
            acf(predicted[ord] - meanRT[ord], lag=500)
            abline(h=0, b=1, col="red", lty=2)
        } else {
            plot(meanRT, predicted, main=sprintf("2 piece linear (%i)", exp), ylab="predicted")
            abline(a=0, b=1, col="red")
        }
##        abline(h=pars[4], col="dark green")
        count <- count + 1
        print(floor(count / 2) + 1)
        par(mfg=c(count %% 2 + 1, floor(count / 2) + 1))
        
        if(residual) {
            ord <- order(meanRT)
            plot(meanRT[ord], mus[ord] - meanRT[ord], main=sprintf("Linear (%i)", exp), ylab="raw", type="l")
            abline(h=0, col="red", lty=2)
        } else if (acf) {
            plot(0, 0)
        } else {
            plot(meanRT, mus, main=sprintf("Linear (%i)", exp), ylab="raw")
            coef <- lm(mus ~ meanRT)$coefficients
            abline(a=coef[1], b=coef[2], col="red")
        }
        count <- count + 1
    }

    if(!is.null(pdf_file))
        dev.off()
}

experiments_list <- split(unique(exp_id), cut(unique(exp_id), breaks=length(unique(exp_id))/3))
count <- 1
for(vec in experiments_list[1:20]) {
    ##checkExperiment(vec, pdf=sprintf("figs/test_%i.pdf", count), acf=TRUE, residual=FALSE)
    checkExperiment(vec, pdf="~/Desktop/tmp.pdf", residual=TRUE)
    count <- count + 1
}


peptideCounts <- sort(table(peptide_id), decreasing=TRUE)
dev.off()
quartz()
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

tib <- as.tibble(cbind(rt=retention_times, exp=experiment_id, pep=peptide_id))
tib %>% filter(pep==12565) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100)
tib %>% filter(pep==15133) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100)
tib %>% filter(pep==9172) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100)
tib %>% filter(pep==13919) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100)
tib %>% filter(pep==8444) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100)

tib %>% filter(pep %in% c(8444, 12565, 15133, 9172, 13919)) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100) + facet_wrap(~ pep)

tib %>% filter(pep %in% c(8444, 12565, 15133, 9172, 13919)) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100) + facet_wrap(~ pep)

peps <- names(sample(sample(sort(table(tib$pep), decreasing=TRUE)[1:1000]), 12))
tib %>% filter(pep %in% peps) %>% ggplot(aes(rt, fill=as.factor(exp))) + geom_histogram(bins=100) + facet_wrap(~ pep)

tib %>% filter(pep %in% peps) %>% ggplot(aes(x=exp, y=rt, col=as.factor(pep))) + geom_point()

tmp <- evidence[evidence$PEP > .05 & evidence$PEP < 0.06, c("Peptide ID", "Raw file", "Retention time")]
tmp[, c("Peptide ID", "Raw file")]




pep_exp_counts <- table(peptide_id, experiment_id)
countIndices <- sort(as.numeric(pep_exp_counts), index.return=TRUE, decreasing=TRUE)$ix

## column is experiment
cidx <- floor(countIndices / nrow(pep_exp_counts)) + 1
## row is peptide
ridx <- countIndices %% nrow(pep_exp_counts)

for(i in 1:10) {
    rts <- retention_times[peptide_id == ridx[i] & experiment_id == cidx[i]]
}
## stanRes <- stan(file="fit_RT.stan", data=data)
