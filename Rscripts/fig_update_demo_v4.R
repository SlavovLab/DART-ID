source('Rscripts/lib.R')

## attempt #3 ----------

set.seed(3)

exps <- 1:6
# generated with rnorm(10, 1, 0.1)
#exp_transforms <- c(1.02, 1.03, 0.98, 0.97, 1.01, 0.99, 1.01, 1.02, 0.96, 0.99)
exp_transforms <- rnorm(length(exps), 1, 0.1)
#exp_transforms <- exp_transforms + rnorm(length(exps), 0, 0.0005)
#exp_transforms

e1_beta <- c(-2, 0.6, 1.4, 20)
e2_beta <- c(2, 1.3, 0.7, 25)

mus <- c(6, 17, 25, 40)
musx <- rep(mus, each=2)
two_piece_adjust <- function(x, e) {
  if(x < e[4]) {
    e[1] + (x * e[2])
  } else {
    e[1] + (e[2] * e[4]) + ((x-e[4])*e[3])
  }
}
rts <- as.vector(sapply(mus, function(x) {
  c(two_piece_adjust(x, e1_beta), two_piece_adjust(x, e2_beta))
}))
rts <- rts + rnorm(length(rts), 0, 1)

xmin <- 0
xmax <- 45

plot(musx, rts, pch=16)
segments(x0=0, x1=e1_beta[4],
         y0=e1_beta[1], y1=e1_beta[1] + (e1_beta[2]*e1_beta[4]))
segments(x0=e1_beta[4], x1=xmax,
         y0=e1_beta[1] + (e1_beta[2]*e1_beta[4]), 
         y1=e1_beta[1] + (e1_beta[2]*e1_beta[4]) + ((xmax-e1_beta[4]) * e1_beta[3]))
segments(x0=0, x1=e2_beta[4],
         y0=e2_beta[1], y1=e2_beta[1] + (e2_beta[2]*e2_beta[4]))
segments(x0=e2_beta[4], x1=xmax,
         y0=e2_beta[1] + (e2_beta[2]*e2_beta[4]), 
         y1=e2_beta[1] + (e2_beta[2]*e2_beta[4]) + ((xmax-e2_beta[4]) * e2_beta[3]))