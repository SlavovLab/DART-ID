source('process.desc.R')
source('parse.ev.adj.R')
#source('adjust.pep.ali.R')
source('adjust.pep.expcentric.R')

load.ev <- function(ev, par.file='dat/params.Fit2.RData', include.REV=FALSE, include.CON=FALSE) {
  load(par.file)
  # remove abnormal LC experiments
  # load experiments from correlation testing in similar.lc.R
  exps.lc <- unlist(read_csv('dat/exps.corr.txt')[,2])
  names(exps.lc) <- NULL
  
  ev <- ev %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>% # Only use Elite experiments
    filter(`Raw file` %in% exps.lc) # Remove abnormal LC experiments
  
  ## Filter of PEP < .05
  ev.f <- ev %>% filter(PEP < 0.05) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) # Only use Elite experiments

  # Remove Reverse matches
  if(!include.REV) {
    ev.f <- ev.f %>% filter(!grepl('REV*', `Leading razor protein`)) 
  }  
  # Remove Contaminants
  if(!include.CON) {
    ev.f <- ev.f %>% filter(!grepl('CON*',`Leading razor protein`))
  }
  
  ev.f <- ev.f %>% 
    filter(PEP < 0.05) %>%
    mutate(Protein=`Leading razor protein`) %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP", "Protein") %>%
    mutate(exp_id=`Raw file`) %>%  # new column - exp_id = numeric version of experiment file
    mutate_at("exp_id", funs(as.numeric(as.factor(.)))) %>%
    mutate(`Stan ID`=`Peptide ID`) %>%
    mutate_at("Stan ID", funs(as.numeric(as.factor(.))))
  ev.f <<- ev.f
  
  experiment_factors <<- as.factor(ev.f$`Raw file`)
  num_exps <<- length(unique(ev.f[["exp_id"]]))
  
  ## "true peptide id" matches peptide id in evidence file
  ## "peptide id" is index in 1:num_peptides for stan
  raw_peptide_id <<- ev.f[["Peptide ID"]]
  pep.id.list <<- unique(raw_peptide_id)
  num_peptides <<- length(unique(ev.f$`Stan ID`))
  
  # parse linear regression params from STAN output
  beta0 <<- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
  beta1 <<- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
  beta2 <<- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
  
  # sigma params from Fit2
  split.point <<- pars[sprintf('split_point[%i]', seq(1, num_exps))]
  sigma.slope <<- pars[sprintf('sigma_slope[%i]', seq(1, num_exps))]
  sigma.intercept <<- pars[sprintf('sigma_intercept[%i]', seq(1, num_exps))]
  sigma.slope.global <<- pars['sigma_slope_global']
  
  # sigma params from Fit1
  sigmas <<- pars[grep('sigma\\[', names(pars))]
  
  mus <<- pars[grep('mu\\[', names(pars))]
}