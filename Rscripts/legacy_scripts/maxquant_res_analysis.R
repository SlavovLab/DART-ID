library(tidyverse)
library(ggridges)
library(RColorBrewer)
source('Rscripts/lib.R')

## load data from file --------------------------------------------------

# generated 20180712
#load('dat/error_df.RData')
load('dat/error_df_20180716.RData')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ev_updated.txt')

ev.f <- ev %>%
  filter(!is.na(pep_new)) %>%
  filter(!grepl('SQC9', `Raw file`)) %>%
  filter(PEP < 0.01)

# coefs <- lm(ev.f$`Retention time` ~ ev.f$`Calibrated retention time`)$coefficients
# ev.f$RT_calibrated <- ((ev.f$`Calibrated retention time` * coefs[2]) + coefs[1])
# 
# plot(ev.f$`Retention time`[1:1e4], ev.f$RT_calibrated[1:1e4])
# abline(a=0, b=1, col='red')

exps <- sort(unique(ev.f$`Raw file`))
ev.mq <- data.frame()
for(exp in exps) {
  ev.a <- ev.f %>% filter(PEP < 0.01) %>% filter(`Raw file`==exp) %>%
    filter(!is.na(`Calibrated retention time`))
  if(nrow(ev.a) == 0) return(c(0))
  coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
  #ev.a$`Retention time` - ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.mq <- rbind(ev.mq, ev.a)
}

ev.mq <- ev.mq %>% arrange(id)
ev.f$RT_calibrated <- ev.mq$RT_calibrated


# plot outliers from DART -------------------------------------------------

ev.f <- ev.f %>%
  mutate(dart_error=`Retention time`-muij,
         mbr_error=`Retention time`-RT_calibrated)

ev.do <- ev.f %>% filter(abs(dart_error) > 5) %>%
  dplyr::select(c('Sequence', 'Raw file', 'PEP', 'Retention time', 'RT_calibrated', 'muij', 'id'))

plot(ev.f$muij[1:1e4], ev.f$RT_calibrated[1:1e4], col=rgb(0,0,0,0.1), pch=16)
abline(a=0, b=1, col='red')

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
k <- 50
contour_cols <- viridis::viridis(k, alpha=0.3)
dens <- get_density(ev.f$dart_error[1:1e4], ev.f$mbr_error[1:1e4], k)
plot(ev.f$dart_error[1:1e4], ev.f$mbr_error[1:1e4], 
     pch=16, cex=0.5,
     xlim=c(-5, 5), xlab='DART', ylab='MaxQuant', main='Residual RT',
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))])

x <- abs(ev.f$dart_error[1:1e4])
y <- abs(ev.f$mbr_error[1:1e4])
dens <- get_density(x, y, k)
plot(x, y,  pch=16,
     xlim=c(0, 5), xlab='DART', ylab='MaxQuant', main='abs(Residual RT)',
     col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))])
