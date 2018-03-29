library(tidyverse)
library(rmutil)
library(viridis)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.RTL.txt')

ev <- ev %>%
  rename(`Sequence ID`=`Peptide ID`) %>%
  rename(`Peptide ID`=`Mod. peptide ID`) # alias the modified peptide ID as the peptide ID

## Filter of PEP < pep_thresh and remove REV and CON and experiments on wrong LC column
pep_thresh <- 0.5

# load exclusion list - common contaminants + keratins
exclude <- read_lines('pd_exclude.txt')

ev.f <- ev %>% 
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
toRemove <- ev.f %>% 
  group_by(`Peptide ID`) %>% 
  summarise(num.exps=length(unique(`exp_id`))) %>% 
  filter(num.exps <= 10) %>%
  pull(`Peptide ID`)
ev.f <- ev.f %>% 
  filter(!(`Peptide ID` %in% toRemove))

exp_id <- as.numeric(as.factor(ev.f$`Raw file`))
exp_names <- levels(as.factor(ev.f$`Raw file`))
num_experiments = length(unique(exp_id))

raw_peptide_id <- ev.f[["Peptide ID"]]
stan_peptide_id <- as.numeric(as.factor(ev.f[["Peptide ID"]]))

retention_times <- ev.f$`Retention time`
pep <- ev.f$PEP

# Initialize canonical retention time priors for each peptide
muInit <- sapply(unique(stan_peptide_id), function(pid) {
  # mean of peptide retention times, weighted by PEP
  weights <- ((1 - pep[stan_peptide_id == pid]) - (1 - pep_thresh)) / pep_thresh
  sum(retention_times[stan_peptide_id==pid] * weights) / sum(weights)
  #sum(retention_times[stan_peptide_id==pid] * weights) / sum(weights) + rnorm(1, 0, 10)
})

# initialize priors for the segmented linear regression
# first element of vector is beta_0, or the intercept
# second element is beta_1 and beta_2, the slopes of the two segments
beta_init <- rbind(rep(10, num_experiments), rep(1, num_experiments))

for(i in 1:10) {
  lm_coefs <- sapply(sort(unique(exp_id)), function(id) {
    rt_cur <- retention_times[exp_id==id]
    mu_cur <- muInit[match(stan_peptide_id[exp_id == id], unique(stan_peptide_id))]
    pep_cur <- pep[exp_id==id]
    lm_cur <- lm(rt_cur ~ mu_cur, weights = 1-pep_cur)
    lm_cur$coefficients
  })
  
  # store linear regression parameters
  beta_init[1, ] <- lm_coefs[1,]
  beta_init[2, ] <- lm_coefs[2,]
  
  # store this set of priors before iterating again
  muInitPrev <- muInit
  
  # calculate new set of canonical RTs based on linear regression params
  mu_pred <- (retention_times - beta_init[1, exp_id]) / beta_init[2, exp_id] 
  # make sure new canonical RTs are within same range as distorted RTs
  mu_pred[mu_pred <= 0] <- min(retention_times)
  mu_pred[mu_pred >= max(retention_times)] <- max(retention_times)
  
  muPrev <- muInit
  
  # new set of priors for canonical RTs based on weighted combination of
  # this set of predicted canonical RTs
  muInit <- sapply(unique(stan_peptide_id), function(pid) {
    weights <- ((1-pep[stan_peptide_id==pid]) - (1-pep_thresh))/pep_thresh
    sum(mu_pred[stan_peptide_id==pid] * weights) / sum(weights)
  })
  
  print(sum(muPrev - muInit)^2 / length(muInit))
}



for(id in 1:max(exp_id)) {
  
  pdf(file=paste('exp_fits_linear/', id, '.pdf', sep=''), width=5, height=5)
  
  plot(muInit[match(stan_peptide_id[exp_id == id], unique(stan_peptide_id))], 
       retention_times[exp_id==id],
       xlab='mu', ylab='rt', 
       main=paste('exp: ', id, ' | n=', sum(exp_id==id), '\n',
                  'b0=',format(beta_init[1,id],digits=3), 
                  ' | b1=', format(beta_init[2,id], digits=3), sep=''))
  abline(a=beta_init[1,id], b=beta_init[2,id], col='red')
  
  dev.off()
}



###


load('dat/params.Fit.linear.RData')
lpars <- pars

load('dat/params.Fit3c.RTL.RData')
spars <- pars

rm(pars)

lmuij_fit <- lpars[grep("muij", names(lpars))]
smuij_fit <- spars[grep("muij", names(spars))]
lmus <- lpars[grep("mu\\[", names(lpars))]
smus <- spars[grep("mu\\[", names(spars))]
lsigma_ij <- lpars[grep("sigma_ij", names(lpars))]
ssigma_ij <- spars[grep("sigma_ij", names(spars))]

pep_exp_all <- paste(stan_peptide_id, exp_id, sep=" - ")
pep_exp_pairs <- unique(pep_exp_all)

muij_map <- match(paste(stan_peptide_id, exp_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

pep_col_code <- cut(pep, breaks=10)

for(exp in 1:num_experiments) {
  print(exp)
  
  sbetas <- c(spars[paste0("beta_0[", exp, "]")], 
              spars[paste0("beta_1[", exp, "]")], 
              spars[paste0("beta_2[", exp, "]")],
              spars[paste0("split_point[", exp, "]")])
  lbetas <- c(lpars[paste0("beta_0[", exp, "]")], 
              lpars[paste0("beta_1[", exp, "]")])
  
  ## get experiment indices
  exp_indices <- muij_to_exp==exp
  
  observed <- retention_times[exp_indices[muij_map]]
  lmu <- lmus[(muij_to_pep[muij_map])[exp_indices[muij_map]]]
  smu <- smus[(muij_to_pep[muij_map])[exp_indices[muij_map]]]
  obs_code <- pep_col_code[exp_indices[muij_map]]
  
  pdf(file=paste0('linear_v_segmented_figs/exp', exp, '.pdf'), width=8, height=10)
  par(mfrow=c(3,2), oma=c(0, 0, 2, 0))
  
  plot(lmu, observed, pch=19, cex=0.2, 
       xlab='Canonical RT', ylab='Observed RT',
       main='Linear Fit')
  abline(a=lbetas[1], b=lbetas[2], col='green')
  
  plot(smu, observed, pch=19, cex=0.2,
       xlab='Canonical RT', ylab='Observed RT',
       main='Segmented Fit')
  abline(v=sbetas[4], col="blue", lty=2)
  segments(x0=0, y0=sbetas[1], x1=sbetas[4], y1=sbetas[1]+sbetas[2]*sbetas[4], col="green", lwd=1.5)
  segments(x0=sbetas[4], y0=sbetas[1]+sbetas[2]*sbetas[4], x1=400, y1=sbetas[1] + sbetas[2]*sbetas[4] +sbetas[3]*(400-sbetas[4]), col="red", lwd=1.5)
  
  lpredicted <- (lmuij_fit[muij_map])[exp_indices[muij_map]]
  spredicted <- (smuij_fit[muij_map])[exp_indices[muij_map]]
  lresiduals <- observed - lpredicted
  sresiduals <- observed - spredicted
  lpredicted_sd <- (lsigma_ij[muij_map])[(muij_to_exp==exp)[muij_map]]
  spredicted_sd <- (ssigma_ij[muij_map])[(muij_to_exp==exp)[muij_map]]
  
  cols <- rev(inferno(n=10))
  
  plot(lpredicted, lresiduals, pch=19, cex=0.2, col=cols[obs_code],
       xlab='Predicted RT', ylab='Residual RT')
  lines(lpredicted[order(lpredicted)],
        sapply(lpredicted_sd, function(s) qlaplace(.025, 0, s))[order(lpredicted)],
        col="red")
  lines(lpredicted[order(lpredicted)],
        sapply(lpredicted_sd, function(s) qlaplace(.975, 0, s))[order(lpredicted)],
        col="red")
  
  plot(spredicted, sresiduals, pch=19, cex=0.2, col=cols[obs_code],
       xlab='Predicted RT', ylab='Residual RT')
  lines(spredicted[order(spredicted)],
        sapply(spredicted_sd, function(s) qlaplace(.025, 0, s))[order(spredicted)],
        col="red")
  lines(spredicted[order(spredicted)],
        sapply(spredicted_sd, function(s) qlaplace(.975, 0, s))[order(spredicted)],
        col="red")
  
  plot(density(lresiduals, n=1e4), xlim=c(-15, 15), ylim=c(0, 0.5), col='red',
       main='Residuals')
  lines(density(sresiduals, n=1e4), col='blue')
  legend(5, 0.2, c('Linear', 'Segmented'), col=c('red', 'blue'), 
         lty = c(1, 1), pch = c(NA, NA), cex=0.8,
         merge=T)
  
  plot(seq(-15,15,by=0.01), ecdf(lresiduals)(seq(-15,15,by=0.01)), type='l', col='red',
       xlab='Residuals', ylab='Density', main='ECDF')
  lines(seq(-15,15,by=0.01), ecdf(sresiduals)(seq(-15,15,by=0.01)), type='l', col='blue')
  legend(5, 0.4, c('Linear', 'Segmented'), col=c('red', 'blue'), 
         lty = c(1, 1), pch = c(NA, NA), cex=0.8,
         merge=T)
  
  mtext(paste0(exp_names[exp]), outer = TRUE, cex = 1.5)
  
  dev.off()
}


# overall residuals
lresiduals <- ev.f$`Retention time`-lmuij_fit[muij_map]
sresiduals <- ev.f$`Retention time`-smuij_fit[muij_map]

plot(density(sresiduals, n=1e4), col='blue', xlim=c(-20, 20), ylim=c(0, 0.4))
lines(density(lresiduals, n=1e4), col='red')

plot(seq(-15,15,by=0.01), ecdf(lresiduals)(seq(-15,15,by=0.01)), type='l', col='red',
     xlab='Residuals', ylab='Density', main='ECDF')
lines(seq(-15,15,by=0.01), ecdf(sresiduals)(seq(-15,15,by=0.01)), type='l', col='blue')
legend(5, 0.4, c('Linear', 'Segmented'), col=c('red', 'blue'), 
       lty = c(1, 1), pch = c(NA, NA), cex=0.8,
       merge=T)
