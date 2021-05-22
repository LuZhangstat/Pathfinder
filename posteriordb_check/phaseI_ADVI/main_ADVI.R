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
library(cmdstanr)
source("../utils/sim.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# ADVI settings #
#ADVI_iter = 40000 
#tol_rel_obj = 0.001
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234
M = 20
seed_list = 1:M

load(file = "../results/lp_posteriordb_LBFGS_h10.RData")

# preallocate results #
ADVI_meanfield_draw <- list() 
ADVI_meanfield_center <- list() 
iter_fit_mean <- c()
calls_lp_mean <- c()
calls_gr_mean <- c()


eta_sequence <- c(100, 10, 1, 0.1, 0.01)
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
  
  # set up initials
  data = get_data(po)
  posterior <- to_posterior(model, data)
  init = lapply(initial_ls[[i]], f <- function(x){constrain_pars(posterior, x)})
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  ###  run ADVI through cmdstanr  ###
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ## ADVI test ##
  for(j in 1:M){
    cat(j, "\t")
    ADVI_work <- TRUE
    fit_ADVI <- mod$variational(data = data,
                                seed = seed_list[j],
                                refresh = 0,
                                init = list(init[[j]]),
                                sig_figs = 16,
                                algorithm = "meanfield",
                                save_latent_dynamics = TRUE)#,
                                #iter = ADVI_iter,
                                #tol_rel_obj = tol_rel_obj)
    tryCatch(
      fit_ADVI$lp(),
      error = function(e) { ADVI_work <<- FALSE})
    if(ADVI_work){
      print(summary(fit_ADVI$lp()))
      unconstrained_draws <- unconstrain_cmdstan_vb(model, data, fit_ADVI)
      ADVI_meanfield_draw[[i]] <- unconstrained_draws # record the last sample of phase I
      ADVI_meanfield_center[[i]] <-  colMeans(unconstrained_draws)
      latent_dyn <- read.csv(fit_ADVI$latent_dynamics_files())
      iter_fit <- as.numeric(latent_dyn$X..stan_version_major...2[
        length(latent_dyn$X..stan_version_major...2)-2])
      iter_fit_mean <- c(iter_fit_mean, iter_fit)
      eta_i = which(eta_sequence == fit_ADVI$metadata()$eta)
      calls_lp_mean <- c(calls_lp_mean, iter_fit + 100 * min(eta_i+2, 6))
      calls_gr_mean <- c(calls_gr_mean, iter_fit + 50 * min(eta_i+1, 5))
      break
    }
  }
}

## check fullrank ##
ADVI_fullrank_draw <- list() 
iter_fit_full <- c()
calls_lp_full <- c()
calls_gr_full <- c()
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
  
  # set up initials
  data = get_data(po)
  posterior <- to_posterior(model, data)
  init = lapply(initial_ls[[i]], f <- function(x){constrain_pars(posterior, x)})
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  ###  run ADVI through cmdstanr  ###
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ## ADVI test ##
  for(j in 1:M){
    cat(j, "\t")
    ADVI_work <- TRUE
    fit_ADVI <- mod$variational(data = data,
                                seed = seed_list[j],
                                refresh = 0,
                                init = list(init[[j]]),
                                sig_figs = 16,
                                algorithm = "fullrank",
                                save_latent_dynamics = TRUE)#,
                                #iter = ADVI_iter,
                                #tol_rel_obj = tol_rel_obj)
    tryCatch(
      fit_ADVI$lp(),
      error = function(e) { ADVI_work <<- FALSE})
    if(ADVI_work){
      print(summary(fit_ADVI$lp()))
      unconstrained_draws <- unconstrain_cmdstan_vb(model, data, fit_ADVI)
      ADVI_fullrank_draw[[i]] <- unconstrained_draws # record the last sample of phase I
      latent_dyn <- read.csv(fit_ADVI$latent_dynamics_files())
      iter_fit_full <- c(iter_fit_full, 
                         as.numeric(latent_dyn$X..stan_version_major...2[
                           length(latent_dyn$X..stan_version_major...2)-2]))
      
      latent_dyn <- read.csv(fit_ADVI$latent_dynamics_files())
      iter_fit <- as.numeric(latent_dyn$X..stan_version_major...2[
        length(latent_dyn$X..stan_version_major...2)-2])
      iter_fit_full <- c(iter_fit_full, iter_fit)
      eta_i = which(eta_sequence == fit_ADVI$metadata()$eta)
      calls_lp_full <- c(calls_lp_full, iter_fit + 100 * min(eta_i+2, 6))
      calls_gr_full <- c(calls_gr_full, iter_fit + 50 * min(eta_i+1, 5))
      break
    }
  }
}


save(file = "../results/ADVI_results_h10.RData",
     list = c("ADVI_meanfield_draw", "ADVI_meanfield_center",
              "ADVI_fullrank_draw", "calls_lp_mean", "calls_gr_mean",
              "calls_lp_full", "calls_gr_full"))



### 

fit_ADVI$cmdstan_diagnose()

fit_ADVI <- mod$variational(data = data,
                            seed = seed_list[j],
                            #refresh = 1,
                            init = list(init[[j]]),
                            sig_figs = 16,
                            algorithm = "meanfield",
                            save_latent_dynamics = TRUE,
                            #iter = ADVI_iter,
                            #tol_rel_obj = tol_rel_obj
                            )



