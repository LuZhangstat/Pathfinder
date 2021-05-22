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
M = 5
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234
seed_list = 1:5

load(file = "../results/lp_posteriordb_LBFGS.RData")

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
  if(modelname == "eight_schools-eight_schools_noncentered"){
    INV <- lp_Int_q_posteriordb_8school_noncen(po, alpha)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    INV <- lp_Int_q_posteriordb_gp_pois_regr(po, alpha)
  } else {INV <- lp_Int_q_posteriordb(po, alpha)}
  lp_INV[, i] = INV[1:2]
  lp_mean[i] = INV[3]
  
  data = get_data(po)
  ## ADVI test ##
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  for(j in 1:M){
    ADVI_work <- TRUE
    tryCatch(
      log <- capture.output(
        suppressWarnings(
          phiI_ADVI <- vb(model, data = data, 
                          seed = seed_list[j],
                          algorithm ="meanfield",
                          iter = 40000, 
                          tol_rel_obj = 0.002))), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    tryCatch(
      log <- capture.output(
        ADVI_lps <- lp_recover(model, data = data, 
                               pos_draws = phiI_ADVI@sim$samples)), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    
    if(ADVI_work){
      #phiI_ADVI
      lp_ADVI_INV[, i] = quantile(ADVI_lps, c(alpha/2, 1-alpha/2))
      lp_ADVI_mean[i] <- mean(ADVI_lps)
      cat("ADVI lp INV:", round(lp_ADVI_INV[, i], digits = 3), 
          "ADVI lp_ mean:", round(lp_ADVI_mean[i], digits = 3), "\t",
          "pareto k:", phiI_ADVI@sim$diagnostics$psis$pareto_k  , "\n")}
  }
  
}



phiI_ADVI <- vb(model, data = data, 
                #seed = 2,
                algorithm ="meanfield", iter = 2000, 
                tol_rel_obj = 0.01, 
                grad_samples = 1,
                importance_resampling = TRUE)


phiI_ADVI <- vb(model, data = data, 
                algorithm ="fullrank", iter = 20000, 
                tol_rel_obj = 0.001, 
                grad_samples = 1,
                importance_resampling = TRUE)

phiI_ADVI@sim$diagnostics$psis$pareto_k
ADVI_lps <- lp_recover(model, data = data, 
                       pos_draws = phiI_ADVI@sim$samples)
quantile(ADVI_lps, c(alpha/2, 0.5, 1-alpha/2))
mean(ADVI_lps)
INV

round(quantile(lp_f_draws, c(alpha/2, 1-alpha/2)), digits = 2)

## Check fullrank ##
# preallocate results #
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
lp_ADVI_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_ADVI_mean <- c()

which(model_record ==24)
for(i in 9:length(model_record)){
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
  if(modelname == "eight_schools-eight_schools_noncentered"){
    INV <- lp_Int_q_posteriordb_8school_noncen(po, alpha)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    INV <- lp_Int_q_posteriordb_gp_pois_regr(po, alpha)
  } else {INV <- lp_Int_q_posteriordb(po, alpha)}
  lp_INV[, i] = INV[1:2]
  lp_mean[i] = INV[3]
  
  data = get_data(po)
  ## ADVI test ##
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  for(j in 1:M){
    ADVI_work <- TRUE
    tryCatch(
      log <- capture.output(
        suppressWarnings(
          phiI_ADVI <- vb(model, data = data, 
                          seed = seed_list[j],
                          algorithm ="fullrank",
                          iter = 40000, 
                          tol_rel_obj = 0.002,
                          importance_resampling = TRUE))), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    tryCatch(
      log <- capture.output(
        ADVI_lps <- lp_recover(model, data = data, 
                               pos_draws = phiI_ADVI@sim$samples)), # meanfield fullrank, 
      error = function(e) { ADVI_work <<- FALSE})
    
    if(ADVI_work){
      #phiI_ADVI
      lp_ADVI_INV[, i] = quantile(ADVI_lps, c(alpha/2, 1-alpha/2))
      lp_ADVI_mean[i] <- mean(ADVI_lps)
      cat("ADVI lp INV:", round(lp_ADVI_INV[, i], digits = 3), 
          "ADVI lp_ mean:", round(lp_ADVI_mean[i], digits = 3), "\t",
          "pareto k:", phiI_ADVI@sim$diagnostics$psis$pareto_k  , "\n")}
  }
  
}

# test code #
phiI_ADVI <- vb(model, data = data, tol_rel_obj = 0.01,
                seed = 123, 
                algorithm ="fullrank")
print(phiI_ADVI, c("mu", "tau", "phi", "E_deaths"))
ADVI_lps <- lp_recover(model, data = data, 
                     pos_draws = phiI_ADVI@sim$samples)
Hk = cov(phiI_ADVI@sim$samples[[1]])

quantile(ADVI_lps, c(alpha/2, 1-alpha/2))
mean(ADVI_lps)


ADVI_fs <- f_recover(model, data = data, 
                       pos_draws = phiI_ADVI@sim$samples)

quantile(ADVI_fs, c(alpha/2, 1-alpha/2))
mean(ADVI_fs)
