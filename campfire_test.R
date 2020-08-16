rm(list = ls())

## test with campfire from ben ##
devtools::install_github("bbbales2/campfire")
library(cmdstanr)
library(posterior)
library(campfire)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


