
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
  pos_draws <- reference_posterior_draws(po, type = "draws")
  lp_recovers <- lp_recover(model, data, pos_draws)
  quantile(lp_recovers, probs = c(alpha / 2.0, 1.0 - alpha / 2.0))
  
}

lp_recover <- function(model, data, pos_draws){
  
  ###
  #' recover the log density based on the posterior draws
  #' 
  ###
  
  posterior <- to_posterior(model, data)
  fn <- function(theta){
    log_prob(posterior, unconstrain_pars(posterior, theta), 
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
  sapply(phI_sample@sim$samples, f)
  
}

lp_explore <- function(phiI_sample, INV, L, M){
  
  # get number of iterations for lp__ to reach INV
  lp_phI <- ls_lp_phI(phiI_sample, L)
  lp_in_INV <- (lp_phI > INV[1] & lp_phI < INV[2])
  n_iters <- apply(lp_in_INV, 2, f <- function(x){which(x == TRUE)[1]})
  
  # get corresponding sum of leapfrog numbers
  sampler_params <- get_sampler_params(phI_sample)
  n_sum_leapfrog <- 
    sapply(1:M, g <- function(x){
      sum(sampler_params[[x]][1:n_iters[x], "n_leapfrog__"])})
  
  return(list(n_iters = n_iters, n_sum_leapfrog = n_sum_leapfrog))
}
