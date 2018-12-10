library(tidyverse)
library(pracma)
library(ggridges)
library(RColorBrewer)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_1000cell_20180607_6/ev_updated.txt")
ev_a <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_set1_20180706_2/ev_updated.txt')
ev_b <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_set2_20180706_1/ev_updated.txt')
#ev_all <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_allsets_20180706_1/ev_updated.txt')
#ev_a <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_set_a_20180805/ev_updated.txt')
#ev_b <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_set_b_20180805/ev_updated.txt')
#ev_all <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_all_20180805/ev_updated.txt')

#ev_a <- read_tsv('/gd/bayesian_RT/Alignments/SQC_validation_a_20180807/ev_updated.txt')
#ev_b <- read_tsv('/gd/bayesian_RT/Alignments/SQC_validation_b_20180807/ev_updated.txt')
#ev_all <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_all_20180805/ev_updated.txt')


# filter and intersect ----------------------------------------------------

conf_thresh <- 1e-2

# find peptides in ev_a which are boosted but have no PSMs with PEP < conf_thresh
peps_a <- ev_a %>%
  filter(!is.na(pep_new) & !is.na(PEP)) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh),
            new=sum(PEP > conf_thresh & pep_new < conf_thresh)) %>%
  filter(spec == 0) %>%
  filter(new > 0) %>%
  arrange(`Modified sequence`)

# find peptides in ev_b which have confident PSMs
peps_b <- ev_b %>%
  filter(!is.na(PEP)) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh)) %>%
  filter(spec > 0) %>%
  arrange(`Modified sequence`)

common_peps <- intersect(peps_a$`Modified sequence`, peps_b$`Modified sequence`)

# ev_c <- data.frame(
#   `Modified sequence`=common_peps,
#   rt_a=as.numeric(peps_a$rt[peps_a$`Modified sequence` %in% common_peps]),
#   rt_b=as.numeric(peps_b$rt[peps_b$`Modified sequence` %in% common_peps]),
#   rt_sd_a=as.numeric(peps_a$rt_sd[peps_a$`Modified sequence` %in% common_peps]),
#   rt_sd_b=as.numeric(peps_b$rt_sd[peps_b$`Modified sequence` %in% common_peps])
# )

ev_c <- ev_a %>% 
  filter(`Modified sequence` %in% common_peps) %>%
  dplyr::select(c('Modified sequence', 'Retention time'))

ev_c$rt_conf <- peps_b$rt[match(ev_c$`Modified sequence`, peps_b$`Modified sequence`)]

# confident observations
conf_peps_a <- unique(ev_a %>% 
  filter(grepl('FP18[A-E]', `Raw file`)) %>%
  filter(PEP < 1e-3) %>%
  pull(`Modified sequence`))
conf_peps_b <- unique(ev_b %>% 
  filter(!grepl('FP18[A-E]', `Raw file`)) %>%
  filter(PEP < 1e-3) %>%
  pull(`Modified sequence`))
conf_peps <- intersect(conf_peps_a, conf_peps_b)

ev_fa <- ev_a %>% 
  filter(grepl('FP18[A-E]', `Raw file`) & PEP < 1e-3) %>%
  filter(`Modified sequence` %in% conf_peps) %>%
  arrange(`Modified sequence`)
ev_fb <- ev_b %>% 
  filter(!grepl('FP18[A-E]', `Raw file`) & PEP < 1e-3) %>%
  filter(`Modified sequence` %in% conf_peps) %>%
  group_by(`Modified sequence`) %>% summarise(rt=median(`Retention time`))

conf_x <- ev_fa$`Retention time`
conf_y <- ev_fb$rt[match(ev_fa$`Modified sequence`, ev_fb$`Modified sequence`)]

# scramble
set.seed(1)
random_x <- ev_fa$`Retention time`
random_y <- ev_fb$rt[sample(match(ev_fa$`Modified sequence`, ev_fb$`Modified sequence`),
                            size=nrow(ev_fa))]

## ----------

pdf(file='manuscript/Figs/rt_validation_v8.pdf', width=3.5, height=2)

#layout(rbind(c(1,1),
#             c(2,3)))

par(pty='s', 
    mar=c(1.75,1.5,0.25,6.5), 
    cex.axis=0.85)

#plot(ev_c$rt_b, ev_c$rt_a, pch=16, cex=0.75,
plot(0, 0, type='n',
     xlim=c(14, 50.5), ylim=c(14, 50.5),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)

#abline(a=0, b=1, col='red')

points(ev_c$rt_conf, ev_c$`Retention time`, pch=16, cex=0.25,
       col=rgb(1,0,0,0.5))

# reduce plotting of conf ones
set.seed(2)
conf_n <- sample.int(length(conf_x), size=300)

points(conf_x[conf_n], conf_y[conf_n], pch=16, col=rgb(0,0,1,0.5), cex=0.25)
points(random_x, random_y, pch=16, col=rgb(0,0,0,0.1), cex=0.25)

axis(1, at=seq(10, 50, by=10), tck=-0.02, mgp=c(0, 0.01, 0))
axis(2, at=seq(10, 50, by=10), tck=-0.02, mgp=c(0, 0.3, 0), las=1)

mtext(expression('Median '*RT[B]*' (min)'), 
      side=1, line=0.95, cex=0.8)
mtext(expression(RT[A]*' (min)'), 
      side=2, line=1.25, cex=0.8)

#'A^1: Boosted peptides\nwithout confident\nspectral PSMs'
par(lheight=0.9)
legend(x=52, y=53, 
       c(expression('Subset '*bold(a[1])*':'),
         'Peptides\nwith confident\nspectra', 
         expression('Subset '*bold(a[2])*':'),
         'Upgraded peptides\nwithout confident\nspectra', 
         'Decoy PSMs'),
       pch=16, pt.cex=1, col=c('blue', NA, 'red', NA, rgb(0,0,0,0.4)),
       xjust=0, yjust=1, bty='n', cex=0.85, y.intersp=c(1, 2.1, 2, 2.2, 2.1), 
       xpd=T)

dev.off()


# boxplot -----------------------------------------------------------------

boxs <- list(
  null=log10(abs(random_y-random_x)),
  #conf=rnorm(1e3, 0, 1),
  new=log10(abs(ev_c$`Retention time`-ev_c$rt_conf)),
  conf=log10(abs(conf_y-conf_x))
)

pdf(file='manuscript/Figs/rt_validation_boxplot_v3.pdf', width=3.5, height=2)
par(mar=c(2.5,4,0.5,1), cex.axis=0.85,
    oma=c(0, 0, 1, 0))

boxplot(boxs, horizontal=T, xlab=NA, ylab=NA,
        col=c(rgb(1,1,1), rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
        #col=c(av[4], av[2], av[1]),
        xaxt='n', yaxt='n', ylim=c(-3.25, 1.6),
        outpch=4, outcol=rgb(0,0,0,0.1), outcex=0.75)
axis(1, at=seq(-16, 10, by=1), tck=-0.02, mgp=c(0, 0.1, 0))
axis(2, at=1:3, 
     labels=c('Decoy', expression('Subset '*A[2]), expression('Subset '*A[1])), 
     las=1, tck=-0.02, mgp=c(0, 0.5, 0))
#mtext('log10 | Residual RT |  (min)', side=1, cex=0.8, line=1.25)
mtext(expression('log'[10]*'  | '*RT[A]-RT[B]*' |'), side=1, cex=0.8, line=1.25)

dev.off()
