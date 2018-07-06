library(tidyverse)
source('Rscripts/lib.R')

## load data --------

source('Rscripts/alignment_comparison.R')

## plot --------

boxplot(dart_error, ylim=c(-3, 3))
