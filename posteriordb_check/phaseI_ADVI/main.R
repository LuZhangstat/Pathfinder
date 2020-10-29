rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
source("../utils/sim.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
M = 3
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234
seed_list = 1:3

# pick models #
model_record <- c(1, 2, 3, 4, 6, 8, 11, 13, 15, 20, 23, 25, 26, 27, 29, 31, 33,
                  34, 36, 37, 40, 41, 43, 47, 51, 55, 61, 94, 95)
N_models <- length(model_record)

# preallocate results #
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
lp_ADVI_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_ADVI_mean <- c()

## check meanfield ##
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  # obtain posterior interval of lp__
  INV <- lp_Int_q_posteriordb(po, alpha)
  lp_INV[, i] = INV[c(1, 2)]
  lp_mean[i] = INV[3] 
  
  dat = get_data(po)
  ## ADVI test ##
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  for(j in 1:M){
    ADVI_work <- TRUE
    tryCatch(
      log <- capture.output(
        suppressWarnings(
          phiI_ADVI <- vb(model, data = dat, 
                          seed = seed_list[j],
                          algorithm ="meanfield"))), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    tryCatch(
      log <- capture.output(
        ADVI_lps <- lp_recover(model, data = dat, 
                               pos_draws = phiI_ADVI@sim$samples)), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    
    if(ADVI_work){
      #phiI_ADVI
      lp_ADVI_INV[, i] = quantile(ADVI_lps, c(alpha/2, 1-alpha/2))
      lp_ADVI_mean[i] <- mean(ADVI_lps)
      cat("ADVI lp INV:", round(lp_ADVI_INV[, i], digits = 3), 
          "ADVI lp_ mean:", round(lp_ADVI_mean[i], digits = 3), "\n")}
  }
  
}


## Check fullrank ##
# preallocate results #
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
lp_ADVI_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_ADVI_mean <- c()

for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  # obtain posterior interval of lp__
  INV <- lp_Int_q_posteriordb(po, alpha)
  lp_INV[, i] = INV[c(1, 2)]
  lp_mean[i] = INV[3] 
  
  dat = get_data(po)
  ## ADVI test ##
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  for(j in 1:M){
    ADVI_work <- TRUE
    tryCatch(
      log <- capture.output(
        suppressWarnings(
          phiI_ADVI <- vb(model, data = dat, 
                          seed = seed_list[j],
                          algorithm ="fullrank"))), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    tryCatch(
      log <- capture.output(
        ADVI_lps <- lp_recover(model, data = dat, 
                               pos_draws = phiI_ADVI@sim$samples)), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    
    if(ADVI_work){
      #phiI_ADVI
      lp_ADVI_INV[, i] = quantile(ADVI_lps, c(alpha/2, 1-alpha/2))
      lp_ADVI_mean[i] <- mean(ADVI_lps)
      cat("ADVI lp INV:", round(lp_ADVI_INV[, i], digits = 3), 
          "ADVI lp_ mean:", round(lp_ADVI_mean[i], digits = 3), "\n")}
  }
  
}
