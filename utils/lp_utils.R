library(tibble)
library(dplyr)
library(tidyr)
### functions for obtaining posterior interval of models in posteriordb
lp_Int_q_posteriordb <- function(pick_model, alpha){
  
  ###
  #' Obtain the lower and upper bound of log density posterior interval of the picked model
  #' alpha: return the bounds of (1.0 - alpha) * 100% posterior interval
  #' 
  ###
  
  sc <- stan_code(pick_model)
  model <- stan_model(model_code = sc)
  data <- get_data(pick_model)
  pos_draws <- reference_posterior_draws(pick_model, type = "draws")
  lp_recovers <- lp_recover(model, data, pos_draws)
  c(quantile(lp_recovers, probs = c(alpha / 2.0, 1.0 - alpha / 2.0)),
    mean(lp_recovers))
}

lp_Int_q_posteriordb_8school_noncen <- function(pick_model, alpha){
  
  ###
  #' Obtain the lower and upper bound of log density posterior interval of the picked model
  #' alpha: return the bounds of (1.0 - alpha) * 100% posterior interval
  #' 
  ###
  
  sc <- stan_code(pick_model)
  model <- stan_model(model_code = sc)
  data <- get_data(pick_model)
  pos_draws <- reference_posterior_draws(pick_model, type = "draws")
  f <- function(pos_sample){
    theta_sam <- sapply(pos_sample, cbind)
    theta_trans =  (theta_sam[, 1:8] - theta_sam[, "mu"]) / theta_sam[, "tau"]
    colnames(theta_trans) = c("theta_trans[1]", "theta_trans[2]", 
                              "theta_trans[3]", "theta_trans[4]", 
                              "theta_trans[5]", "theta_trans[6]", 
                              "theta_trans[7]", "theta_trans[8]")
    theta_out = cbind(theta_trans, theta_sam[, 9:10], theta_sam[, 1:8])
    theta_out = as.data.frame(theta_out)
    return(theta_out)
  }
  pos_draws2 = list()
  for (i in 1:length(pos_draws)){
    pos_draws2[[i]] <- f(pos_draws[[i]])
  }
  lp_recovers <- lp_recover(model, data, pos_draws2)
  c(quantile(lp_recovers, probs = c(alpha / 2.0, 1.0 - alpha / 2.0)),
    mean(lp_recovers))
}


lp_Int_q_posteriordb_gp_pois_regr <- function(pick_model, alpha){
  
  ###
  #' Obtain the lower and upper bound of log density posterior interval of the picked model
  #' alpha: return the bounds of (1.0 - alpha) * 100% posterior interval
  #' 
  ###
  
  sc <- stan_code(pick_model)
  model <- stan_model(model_code = sc)
  data <- get_data(pick_model)
  recover_f_tilde <- function(theta, x, N){
    cov = theta[2]^2 *exp(-0.5 / (theta[1]^2) * as.matrix(dist(x)^2)) + 
      diag(N)*1e-10;
    L_cov = chol(cov);
    f_tilde = forwardsolve(L_cov, theta[3:13], upper.tri = TRUE,
                 transpose = TRUE);
    return(f_tilde);
  }
  pos_draws <- reference_posterior_draws(pick_model, type = "draws")
  f <- function(pos_sample){
    theta_sam <- sapply(pos_sample, cbind)
    theta_trans =  t(apply(theta_sam, 1, recover_f_tilde, data$x, data$N))
    colnames(theta_trans) = c("f_tilde[1]", "f_tilde[2]", 
                              "f_tilde[3]", "f_tilde[4]",
                              "f_tilde[5]", "f_tilde[6]",
                              "f_tilde[7]", "f_tilde[8]",
                              "f_tilde[9]", "f_tilde[10]", "f_tilde[11]")
    theta_out = cbind(theta_sam[, 1:2], theta_trans)
    theta_out = as.data.frame(theta_out)
    return(theta_out)
  }
  pos_draws2 = list()
  for (i in 1:length(pos_draws)){
    pos_draws2[[i]] <- f(pos_draws[[i]])
  }
  lp_recovers <- lp_recover(model, data, pos_draws2)
  c(quantile(lp_recovers, probs = c(alpha / 2.0, 1.0 - alpha / 2.0)),
    mean(lp_recovers))
}


lp_recover <- function(model, data, pos_draws){
  
  ###
  #' recover the log density based on the posterior draws
  #' 
  ###
  
  posterior <- to_posterior(model, data)
  par_template <- get_inits(posterior)[[1]]  # get random inits
  npar <- length(par_template)              # number of pars
  n_inits <- sapply(par_template, length)
  
  fn <- function(theta){
    j = 1
    for(i in 1:npar){
      par_template[[i]] <- theta[j:(j + n_inits[i] - 1)]
      j = j + n_inits[i]
    }
    log_prob(posterior, unconstrain_pars(posterior, par_template), 
             adjust_transform = TRUE, gradient = TRUE)[1]
  }
  lpn <- function(gsd_l) apply(sapply(gsd_l, unlist), 1, fn)
  sapply(pos_draws, lpn)
}

# unconstrain_vb <- function(model, data, vb_obj) {
#   post <- to_posterior(model, data)
#   
#   tmp <- rstan:::create_skeleton(post@model_pars, post@par_dims)
#   skeleton <- tmp[names(tmp) != "lp__"]
#   constrained_draws <- as.matrix(sapply(vb_obj@sim$samples[[1]], unlist))
#   S <- nrow(constrained_draws)
#   unconstrained_draws <- sapply(seq_len(S), FUN = function(s) {
#     constrained_s <- rstan:::rstan_relist(constrained_draws[s, ], skeleton)
#     post@.MISC$stan_fit_instance$unconstrain_pars(constrained_s)
#   })
#   unconstrained_draws <- t(unconstrained_draws)
#   unconstrained_draws
# }


unconstrain_cmdstan_vb <- function(model, data, vb_obj) {
  post <- to_posterior(model, data)
  
  tmp <- rstan:::create_skeleton(post@model_pars, post@par_dims)
  skeleton <- tmp[names(tmp) != "lp__"]
  constrained_draws <- vb_obj$draws()[, c(-1, -2)]
  S <- nrow(constrained_draws)
  unconstrained_draws <- sapply(seq_len(S), FUN = function(s) {
    constrained_s <- rstan:::rstan_relist(constrained_draws[s, ], skeleton)
    post@.MISC$stan_fit_instance$unconstrain_pars(constrained_s)
  })
  unconstrained_draws <- t(unconstrained_draws)
  unconstrained_draws
}


unconstrain_draws <- function(pos_draws, post) {
  post <- to_posterior(model, data)
  
  tmp <- rstan:::create_skeleton(post@model_pars, post@par_dims)
  skeleton <- tmp[names(tmp) != "lp__"]
  constrained_draws <- as.matrix(sapply(pos_draws, unlist))
  S <- nrow(constrained_draws)
  unconstrained_draws <- sapply(seq_len(S), FUN = function(s) {
    constrained_s <- rstan:::rstan_relist(constrained_draws[s, ], skeleton)
    post@.MISC$stan_fit_instance$unconstrain_pars(constrained_s)
  })
  unconstrained_draws <- t(unconstrained_draws)
  unconstrained_draws
}

unconstrain_cmd_draws <- function(last_pI_draws, post) {

  tmp <- rstan:::create_skeleton(post@model_pars, post@par_dims)
  skeleton <- tmp[names(tmp) != "lp__"]
  #constrained_draws <- as.matrix(sapply(pos_draws, unlist))
  constrained_draws <- last_pI_draws[, -1]
  S <- nrow(constrained_draws)
  unconstrained_draws <- sapply(seq_len(S), FUN = function(s) {
    constrained_s <- rstan:::rstan_relist(constrained_draws[s, ], skeleton)
    post@.MISC$stan_fit_instance$unconstrain_pars(constrained_s)
  })
  unconstrained_draws <- t(unconstrained_draws)
  unconstrained_draws
}


modify_8school_noncen <- function(pick_model){
  
  ###
  #' Obtain posterior sample of the picked model
  #' alpha: return the bounds of (1.0 - alpha) * 100% posterior interval
  #' 
  ###
  
  sc <- stan_code(pick_model)
  model <- stan_model(model_code = sc)
  data <- get_data(pick_model)
  pos_draws <- reference_posterior_draws(pick_model, type = "draws")
  f <- function(pos_sample){
    theta_sam <- sapply(pos_sample, cbind)
    theta_trans =  (theta_sam[, 1:8] - theta_sam[, "mu"]) / theta_sam[, "tau"]
    colnames(theta_trans) = c("theta_trans[1]", "theta_trans[2]", 
                              "theta_trans[3]", "theta_trans[4]", 
                              "theta_trans[5]", "theta_trans[6]", 
                              "theta_trans[7]", "theta_trans[8]")
    theta_out = cbind(theta_trans, theta_sam[, 9:10], theta_sam[, 1:8])
    theta_out = as.data.frame(theta_out)
    return(theta_out)
  }
  pos_draws2 = list()
  for (i in 1:length(pos_draws)){
    pos_draws2[[i]] <- f(pos_draws[[i]])
  }
  return(pos_draws2)
}

modify_8school_noncen2 <- function(theta_sam){
  
  theta_trans =  (theta_sam[, 1:8] - theta_sam[, 9]) / theta_sam[, 10]
  theta_out = cbind(theta_trans, theta_sam[, 9:10])
  colnames(theta_out) = c("theta[1]", "theta[2]", "theta[3]", "theta[4]",
                          "theta[5]", "theta[6]", "theta[7]", "theta[8]",
                          "mu", "tau")
  #theta_out = as.data.frame(theta_out)
 
  return(theta_out)
}

modify_draws_gp_pois_regr <- function(pick_model){
  
  ###
  #' Obtain posterior samples for the picked model
  #' alpha: return the bounds of (1.0 - alpha) * 100% posterior interval
  #' 
  ###
  
  sc <- stan_code(pick_model)
  model <- stan_model(model_code = sc)
  data <- get_data(pick_model)
  recover_f_tilde <- function(theta, x, N){
    cov = theta[2]^2 *exp(-0.5 / (theta[1]^2) * as.matrix(dist(x)^2)) + 
      diag(N)*1e-10;
    L_cov = chol(cov);
    f_tilde = forwardsolve(L_cov, theta[3:13], upper.tri = TRUE,
                           transpose = TRUE);
    return(f_tilde);
  }
  pos_draws <- reference_posterior_draws(pick_model, type = "draws")
  f <- function(pos_sample){
    theta_sam <- sapply(pos_sample, cbind)
    theta_trans =  t(apply(theta_sam, 1, recover_f_tilde, data$x, data$N))
    colnames(theta_trans) = c("f_tilde[1]", "f_tilde[2]", 
                              "f_tilde[3]", "f_tilde[4]",
                              "f_tilde[5]", "f_tilde[6]",
                              "f_tilde[7]", "f_tilde[8]",
                              "f_tilde[9]", "f_tilde[10]", "f_tilde[11]")
    theta_out = cbind(theta_sam[, 1:2], theta_trans)
    theta_out = as.data.frame(theta_out)
    return(theta_out)
  }
  pos_draws2 = list()
  for (i in 1:length(pos_draws)){
    pos_draws2[[i]] <- f(pos_draws[[i]])
  }
  return(pos_draws2)
}

modify_gp_pois_regr2 <- function(theta_sam){
  
  recover_f_tilde <- function(theta, x, N){
    cov = theta[2]^2 *exp(-0.5 / (theta[1]^2) * as.matrix(dist(x)^2)) + 
      diag(N)*1e-10;
    L_cov = chol(cov);
    f_tilde = forwardsolve(L_cov, theta[3:13], upper.tri = TRUE,
                           transpose = TRUE);
    return(f_tilde);
  }
  po <- posterior("gp_pois_regr-gp_pois_regr", pdb = pd)
  data <- get_data(po)
  
  theta_trans =  t(apply(theta_sam, 1, recover_f_tilde, data$x, data$N))
  theta_out = cbind(theta_sam[, 1:2], theta_trans)
  colnames(theta_out) = c("rho", "alpha", "f_tilde[1]", "f_tilde[2]", 
                          "f_tilde[3]", "f_tilde[4]",
                          "f_tilde[5]", "f_tilde[6]",
                          "f_tilde[7]", "f_tilde[8]",
                          "f_tilde[9]", "f_tilde[10]", "f_tilde[11]")
  
  return(theta_out)
}

lp_recover_2 <- function(model, data, pos_draws){
  
  ###
  #' recover the log density based on the posterior draws
  #' 
  ###
  
  posterior <- to_posterior(model, data)
  unconstrained_draws <- lapply(pos_draws, unconstrain_draws, posterior)
  lps <- unlist(lapply(unconstrained_draws, f <- function(x){
    apply(x, 1, g <- function(y){
      log_prob(posterior, y, adjust_transform = TRUE, gradient = TRUE)[1]})
  }))
  return(lps)
}


f_recover <- function(model, data, pos_draws){
  
  ###
  #' recover the log density based on the posterior draws
  #' 
  ###
  
  posterior <- to_posterior(model, data)
  par_template <- get_inits(posterior)[[1]]  # get random inits
  npar <- length(par_template)              # number of pars
  n_inits <- sapply(par_template, length)
  
  fn <- function(theta){
    j = 1
    for(i in 1:npar){
      par_template[[i]] <- theta[j:(j + n_inits[i] - 1)]
      j = j + n_inits[i]
    }
    log_prob(posterior, unconstrain_pars(posterior, par_template), 
             adjust_transform = TRUE, gradient = TRUE)[1]
  }
  lpn <- function(gsd_l) apply(sapply(gsd_l, unlist), 1, fn)
  sapply(pos_draws, lpn)
}


ls_lp_phI <- function(phiI_sample, L){
  
  ###
  #' retrieve all posterior samples of lp__ from stan phase I samples 
  #' L the number of iterations for phase I
  ###
  
  f <- function(phI_sam){phI_sam$lp__[1:L]}
  sapply(phiI_sample@sim$samples, f)
  
}

lp_explore <- function(fit, INV, L, M, model, data){
  
  # get number of iterations for lp__ to reach INV
  # lp_phI <- ls_lp_phI(phiI_sample, L)
  # lp_phI <- fit$draws("lp__", inc_warmup = TRUE)[, , 1]
  pos_d <- list()
  fit_d <- fit$draws(inc_warmup = TRUE)
  for(ll in 1:dim(fit_d)[2]){
    pos_d[[ll]] <- as.data.frame(fit_d[, ll, -1])
  }
  lp_phI <- lp_recover(model, data, pos_d)
  lp_in_INV <- (lp_phI >= INV[1] & lp_phI <= INV[2])
  n_iters <- apply(lp_in_INV, 2, f <- function(x){which(x == TRUE)[1]})
  if(sum(is.na(n_iters)) > 0){
    print("not all chains reach the target region")
    n_iters[is.na(n_iters)] <- L
  }
  
  # get corresponding sum of leapfrog numbers
  # sampler_params <- get_sampler_params(phiI_sample)
  # n_sum_leapfrog <- 
  #   sapply(1:M, g <- function(x){
  #     sum(sampler_params[[x]][1:n_iters[x], "n_leapfrog__"])})
  
  fit_sum <- fit$sampler_diagnostics(inc_warmup = TRUE)[, , "n_leapfrog__"]
  n_sum_leapfrog <- sapply(1:M, g <- function(x){
    sum(fit_sum[1:n_iters[x], x, "n_leapfrog__"])
  })
  
  return(list(n_iters = n_iters, n_sum_leapfrog = n_sum_leapfrog))
}

# old code #
# ## function for extract optims and inits ##
# get_init_optim <- function(ind){
#   
#   # remove the failed Pathfinder
#   param_path <- lp_opath[[ind]]$opath
#   check <- sapply(param_path, f <- function(x){
#     work <- TRUE
#     tryCatch(
#       lps <- extract_lps(x), error = function(e) { work <<- FALSE})
#     work
#   })
#   filter_mode <- c(1:length(check))[check]
#   
#   param_path <- param_path[filter_mode]
#   
#   lp_ind = ncol(param_path[[1]]$y)
#   inits <- c()
#   optims <- c()
#   for(l in 1:length(param_path)){
#     inits = rbind(inits, param_path[[l]]$y[1, 1:(lp_ind - 1)])
#     last_ind <- nrow(param_path[[l]]$y)
#     optims = rbind(optims, 
#                    param_path[[l]]$y[last_ind, 1:(lp_ind - 1)])
#   }
#   return(list(inits = inits, optims = optims))
# }
