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
set_cmdstan_path("/home/luzhang/Google Drive/Github/adaptation/cmdstan-lowrank_robust")
library(loo)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# ADVI settings #
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234
M = 100
seed_list = 1:100

load(file = "../results/lp_posteriordb_LBFGS_h6.RData")
#load("../results/ADVI_100.RData")

# preallocate results #
ADVI_meanfield_draw <- list() 
ADVI_meanfield_center <- list() 
ADVI_meanfield_draw_100 <- list()
ADVI_meanfield_lrs_100 <- list()
iter_fit_mean <- array(data = 0, dim = c(M, length(model_record)))
calls_lp_mean <- array(data = 0, dim = c(M, length(model_record)))
calls_gr_mean <- array(data = 0, dim = c(M, length(model_record)))

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
  D <- get_num_upars(posterior)
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  ###  run ADVI through cmdstanr  ###
  # need to run the following line to generate the stan file. no need to run if there is one. 
  # write_stan_file(sc, dir = paste0(getwd(), "/modelcode"),
  #                 basename = paste0(modelname, ".stan"))
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ## 100 ADVI meanfield ##
  ADVI_meanfield_center[[i]] <- array(data = NA, dim = c(D, M))
  ADVI_meanfield_draw[[i]] <- array(data = NA, dim = c(D, M))
  ADVI_meanfield_draw_100[[i]] <- list()
  ADVI_meanfield_lrs_100[[i]] <- list()
  j = 1
  seed = 1
  while(j <= M){
    cat(j, "\t")
    ADVI_work <- TRUE
    fit_ADVI <- mod$variational(data = data,
                                seed = seed,
                                refresh = 0,
                                sig_figs = 18,
                                algorithm = "meanfield",
                                save_latent_dynamics = TRUE,
                                output_samples = 100)
    tryCatch(
      unconstrained_draws <- unconstrain_cmdstan_vb(model, data, fit_ADVI),
      error = function(e) { ADVI_work <<- FALSE})
    if(ADVI_work){
      ADVI_meanfield_draw_100[[i]][[j]] <- unconstrained_draws # record the last sample of phase I
      ADVI_meanfield_center[[i]][, j] <-  colMeans(unconstrained_draws)
      lrs <- c(fit_ADVI$draws("lp__") - fit_ADVI$draws("lp_approx__"))
      ADVI_meanfield_lrs_100[[i]][[j]] <- lrs
      finit_ind <- which(is.finite(lrs))
      pick_ind <- sample(finit_ind, size = 1)
      ADVI_meanfield_draw[[i]][, j] <- unconstrained_draws[pick_ind, ]
      latent_dyn <- read.csv(fit_ADVI$latent_dynamics_files())
      iter_fit <- as.numeric(latent_dyn$X..stan_version_major...2[
        length(latent_dyn$X..stan_version_major...2)-2])
      iter_fit_mean[j, i] <- iter_fit
      eta_i = which(eta_sequence == fit_ADVI$metadata()$eta)
      calls_lp_mean[j, i] <- calls_lp_mean[j, i] + iter_fit + 100 * min(eta_i+2, 6)
      calls_gr_mean[j, i] <- calls_gr_mean[j, i] + iter_fit + 50 * min(eta_i+1, 5)
      j = j + 1
    }else{
      ## if fail, record the estimated number of evaluation of log-density and gradient 
      output_print <- capture.output(fit_ADVI$output())
      #print(fit_ADVI$output())
      L_e <- length(output_print)
      if(output_print[L_e] == "stan::variational::advi::adapt_eta: All proposed step-sizes failed. Your model may be either severely ill-conditioned or misspecified."){
        print("All proposed step-sizes failed")
        calls_lp_mean[j, i] <- calls_lp_mean[j, i] + 100 * 6
        calls_gr_mean[j, i] <- calls_gr_mean[j, i] + 50 * 5
      }else if(output_print[L_e] == "stan::variational::normal_meanfield::calc_grad: The number of dropped evaluations has reached its maximum amount (10). Your model may be either severely ill-conditioned or misspecified." |
               output_print[L_e] == "stan::variational::advi::calc_ELBO: The number of dropped evaluations has reached its maximum amount (100). Your model may be either severely ill-conditioned or misspecified."){
        print("The number of dropped evaluations has reached its maximum amount (10)")
        splited_output = unlist(lapply(output_print, 
                                       f <- function(x){
                                         strsplit(as.character(x),split = " ")}))
        
        eta_ind = which(splited_output == "[eta") + 2 # check line seach times the 
        eta_i = which(eta_sequence == 
                        as.numeric(strsplit(splited_output[[eta_ind]], split = "]")[[1]]))
        calls_lp_mean[j, i] <- calls_lp_mean[j, i] + 100 * min(eta_i+2, 6)
        calls_gr_mean[j, i] <- calls_gr_mean[j, i] + 50 * min(eta_i+1, 5)
        if(output_print[L_e-1] != "  iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes "){
          # started fitting ADVI
          round_iter <- as.numeric(strsplit(
            as.character(output_print[L_e-1]), split = " ")[[1]])
          round_iter <- round_iter[!is.na(round_iter)][1]
          #print(round_iter)
          calls_lp_mean[j, i] <- calls_lp_mean[j, i] + round_iter
          calls_gr_mean[j, i] <- calls_gr_mean[j, i] + round_iter
        }
      }else if(output_print[L_e] == "stan::variational::advi::adapt_eta: Cannot compute ELBO using the initial variational distribution. Your model may be either severely ill-conditioned or misspecified."){
        print("Cannot compute ELBO using the initial variational distribution.")
      }else if(output_print[L_e] == "  iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes "){
        print("stop at the begining of stochastic gradient ascent")
        splited_output = unlist(lapply(output_print, 
                                       f <- function(x){
                                         strsplit(as.character(x),split = " ")}))
        
        eta_ind = which(splited_output == "[eta") + 2 # check line seach times the 
        eta_i = which(eta_sequence == 
                        as.numeric(strsplit(splited_output[[eta_ind]], split = "]")[[1]]))
        calls_lp_mean[j, i] <- calls_lp_mean[j, i] + 100 * min(eta_i+2, 6)
        calls_gr_mean[j, i] <- calls_gr_mean[j, i] + 50 * min(eta_i+1, 5)
      }else{
        print("other problems")
        cat("seed:", seed, "\n")
      }
    }
    seed = seed + 1
  }
}

## check fullrank ##
ADVI_fullrank_draw <- list() 
ADVI_fullrank_draw_100 <- list() 
ADVI_fullrank_lrs_100 <- list()
iter_fit_full <- array(data = 0, dim = c(M, length(model_record)))
calls_lp_full <- array(data = 0, dim = c(M, length(model_record)))
calls_gr_full <- array(data = 0, dim = c(M, length(model_record)))
for(i in 21:length(model_record)){
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
  D <- get_num_upars(posterior)
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  ###  run ADVI through cmdstanr  ###
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ## 100 ADVI meanfield ##
  ADVI_fullrank_draw[[i]] <- array(data = NA, dim = c(D, M))
  ADVI_fullrank_draw_100[[i]] <- list()
  ADVI_fullrank_lrs_100[[i]] <- list()
  j = 1
  seed = 1
  while(j <= M){
    cat(j, "\t")
    ADVI_work <- TRUE
    fit_ADVI <- mod$variational(data = data,
                                seed = seed,
                                refresh = 0,
                                sig_figs = 18,
                                algorithm = "fullrank",
                                save_latent_dynamics = TRUE,
                                output_samples = 100)
    tryCatch(
      unconstrained_draws <- unconstrain_cmdstan_vb(model, data, fit_ADVI),
      error = function(e) { ADVI_work <<- FALSE})
    if(ADVI_work){
      ADVI_fullrank_draw_100[[i]][[j]] <- unconstrained_draws # record the last sample of phase I
      lrs <- c(fit_ADVI$draws("lp__") - fit_ADVI$draws("lp_approx__"))
      ADVI_fullrank_lrs_100[[i]][[j]] <- lrs
      finit_ind <- which(is.finite(lrs))
      pick_ind <- sample(finit_ind, size = 1)
      ADVI_fullrank_draw[[i]][, j] <- unconstrained_draws[pick_ind, ]
      latent_dyn <- read.csv(fit_ADVI$latent_dynamics_files())
      iter_fit <- as.numeric(latent_dyn$X..stan_version_major...2[
        length(latent_dyn$X..stan_version_major...2)-2])
      iter_fit_full[j, i] <- iter_fit
      eta_i = which(eta_sequence == fit_ADVI$metadata()$eta)
      calls_lp_full[j, i] <- calls_lp_full[j, i] + iter_fit + 100 * min(eta_i+2, 6)
      calls_gr_full[j, i] <- calls_gr_full[j, i] + iter_fit + 50 * min(eta_i+1, 5)
      j = j + 1
    }else{
      ## if fail, record the evaluation of log-density and gradient 
      output_print <- capture.output(fit_ADVI$output())
      #print(fit_ADVI$output())
      L_e <- length(output_print)
      if(output_print[L_e] == "stan::variational::advi::adapt_eta: All proposed step-sizes failed. Your model may be either severely ill-conditioned or misspecified."){
        print("All proposed step-sizes failed")
        calls_lp_full[j, i] <- calls_lp_full[j, i] + 100 * 6
        calls_gr_full[j, i] <- calls_gr_full[j, i] + 50 * 5
      }else if(output_print[L_e] == "stan::variational::normal_meanfield::calc_grad: The number of dropped evaluations has reached its maximum amount (10). Your model may be either severely ill-conditioned or misspecified." |
               output_print[L_e] == "stan::variational::advi::calc_ELBO: The number of dropped evaluations has reached its maximum amount (100). Your model may be either severely ill-conditioned or misspecified."){
        print("The number of dropped evaluations has reached its maximum amount (10)")
        splited_output = unlist(lapply(output_print, 
                                       f <- function(x){
                                         strsplit(as.character(x),split = " ")}))
        
        eta_ind = which(splited_output == "[eta") + 2 # check line seach times the 
        eta_i = which(eta_sequence == 
                        as.numeric(strsplit(splited_output[[eta_ind]], split = "]")[[1]]))
        calls_lp_full[j, i] <- calls_lp_full[j, i] + 100 * min(eta_i+2, 6)
        calls_gr_full[j, i] <- calls_gr_full[j, i] + 50 * min(eta_i+1, 5)
        if(output_print[L_e-1] != "  iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes "){
          # started fitting ADVI
          round_iter <- as.numeric(strsplit(
            as.character(output_print[L_e-1]), split = " ")[[1]])
          round_iter <- round_iter[!is.na(round_iter)][1]
          #print(round_iter)
          calls_lp_full[j, i] <- calls_lp_full[j, i] + round_iter
          calls_gr_full[j, i] <- calls_gr_full[j, i] + round_iter
        }
      }else if(output_print[L_e] == "stan::variational::advi::adapt_eta: Cannot compute ELBO using the initial variational distribution. Your model may be either severely ill-conditioned or misspecified."){
        print("Cannot compute ELBO using the initial variational distribution.")
      }else if(output_print[L_e] == "  iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes "){
        print("stop at the begining of stochastic gradient ascent")
        splited_output = unlist(lapply(output_print, 
                                       f <- function(x){
                                         strsplit(as.character(x),split = " ")}))
        
        eta_ind = which(splited_output == "[eta") + 2 # check line seach times the 
        eta_i = which(eta_sequence == 
                        as.numeric(strsplit(splited_output[[eta_ind]], split = "]")[[1]]))
        calls_lp_full[j, i] <- calls_lp_full[j, i] + 100 * min(eta_i+2, 6)
        calls_gr_full[j, i] <- calls_gr_full[j, i] + 50 * min(eta_i+1, 5)
      }else if(output_print[L_e] == "Exception: integrate_ode_rk45:  Failed to integrate to next output time (12) in less than max_num_steps steps (in '/tmp/RtmplArVzk/model-14b74639d7f67.stan', line 32, column 2 to line 35, column 42)"){
        print("integrate_ode_rk45 error")
        
      }else{
        print("other problems")
      }
    }
    seed = seed + 1
  }
}

save(file = "../results/ADVI_100_RVI.RData",
     list = c("ADVI_meanfield_draw", "ADVI_meanfield_center", 
              "ADVI_meanfield_draw_100", "ADVI_meanfield_lrs_100",
              "ADVI_fullrank_draw", "ADVI_fullrank_draw_100",
              "ADVI_fullrank_lrs_100",
              "calls_lp_mean", "calls_gr_mean",
              "calls_lp_full", "calls_gr_full"))



######### old code ###########
# randomly sample 10 ADVI + SIR WOR #
Imp_Resam_WOR_ADVI <- function(ADVI_draws, ADVI_lrs, n_inits, seed = 123){
  
  #' SIR without replacement, index of apporoximating distribution as an auxilary variable
  #' Return n_inits distinct samples 
  
  ## extract samples and log ratios ##
  samples <- do.call("rbind", ADVI_draws)
  lrms <- c(unlist(ADVI_lrs))
  
  ## take off samples with infinite log ratios
  finit_ind <- is.finite(lrms)
  samples <- samples[finit_ind, ]
  lrms <- lrms[finit_ind]
  
  ## compute the importance weight ##
  sample_weights_psis <- suppressWarnings(weights(psis(lrms, r_eff = NA),
                                                  log = FALSE))
  sample_weights_IS <- exp(lrms - max(lrms))/sum(exp(lrms - max(lrms)))
  
  ## Importance resampling ##
  set.seed(seed)
  sample_ind <-sample.int(nrow(samples), size = n_inits, replace = FALSE, 
                          prob = sample_weights_psis)
  
  return(samples[sample_ind, ])
}


MC = 100
load(file = "../results/ADVI_100.RData")
set.seed(123)

ADVI_mf_10_draws <- list()
ADVI_fr_10_draws <- list()

which(model_record == 48)
t_0 <- proc.time()
for(i in 1:49){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  ADVI_mf_10_draws[[i]] <- list()
  ADVI_fr_10_draws[[i]] <- list()
  for(j in 1:MC){
    set.seed(j)
    pick_ind <- sample.int(MC, 10, replace = FALSE)
    ADVI_mf_10_draws[[i]][[j]] <- 
      Imp_Resam_WOR_ADVI(ADVI_meanfield_draw_100[[i]][pick_ind], 
                         ADVI_meanfield_lrs_100[[i]][pick_ind], 
                         n_inits = 10, seed = j)
    ADVI_fr_10_draws[[i]][[j]] <- 
      Imp_Resam_WOR_ADVI(ADVI_fullrank_draw_100[[i]][pick_ind], 
                         ADVI_fullrank_lrs_100[[i]][pick_ind], 
                         n_inits = 10, seed = j)
  }
}
proc.time() - t_0

save(file = "../results/ADVI_10_draws.RData",
     list = c("ADVI_mf_10_draws", "ADVI_fr_10_draws"))





