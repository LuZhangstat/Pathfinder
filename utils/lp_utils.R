
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

lp_explore <- function(phiI_sample, INV, L, M){
  
  # get number of iterations for lp__ to reach INV
  lp_phI <- ls_lp_phI(phiI_sample, L)
  lp_in_INV <- (lp_phI >= INV[1] & lp_phI <= INV[2])
  n_iters <- apply(lp_in_INV, 2, f <- function(x){which(x == TRUE)[1]})
  if(sum(is.na(n_iters)) > 0){
    print("not all chains reach the target region")
    n_iters[is.na(n_iters)] <- L
  }
  
  # get corresponding sum of leapfrog numbers
  sampler_params <- get_sampler_params(phiI_sample)
  n_sum_leapfrog <- 
    sapply(1:M, g <- function(x){
      sum(sampler_params[[x]][1:n_iters[x], "n_leapfrog__"])})
  
  return(list(n_iters = n_iters, n_sum_leapfrog = n_sum_leapfrog))
}
