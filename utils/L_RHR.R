#rm(list = ls())
# L-RHR function 
library(mgcv) # use cholup for one rank Cholesky factor update
#install.packages("rHanso", repos="http://R-Forge.R-project.org")
library(rHanso) 
#help(rHanso)
library(Matrix)

# source("utils/sim_RHR.R")
# source("utils/lp_utils.R")

RHR <- function(fn, 
                gr, 
                x0, 
                sigma0,
                max_iter, 
                abs_tol, 
                accept_tol){
  k = 0; rk = 1; gk = gr(x0); f0 = fn(x0)
  vk = c(sqrt(sum(gk^2)));
  zk = gk / vk[1];
  sigmak = sigma0
  Rk = matrix(sqrt(sigmak), nrow = 1, ncol = 1);
  Zk = matrix(zk, ncol = 1)
  xk = x0; 
  x_list = list(); x_list[[1]] = x0;
  fn_ls = c(); fn_ls[[1]] = f0;
  pk = c(); dk = c(); qk = c()
  Tk = matrix(vk[1], nrow = 1, ncol = 1)
  D = length(gk)
  
  
  # steps and conds of line search
  steps = c(); conds = c();
  
  f_old = f0
  for(k in 1:max_iter){
    dk = forwardsolve(Rk, -vk, upper.tri = TRUE, transpose = TRUE)
    qk = backsolve(Rk, dk, upper.tri = TRUE, transpose = FALSE)
    pk = Zk %*% qk
    pk = pk / sqrt(sum(pk^2)) # different from the paper, for numerical stability
    # Find alphak satisfying the Wolfe conditions (2.2)
    wolfe_cond_res <- linesch_sw(fn = fn,
                                 gr = gr,
                                 xk,
                                 d = pk,
                                 f0 = fn(xk),
                                 grad0 = gr(xk),
                                 c1 = 0, c2 = 0.5,
                                 fvalquit = -Inf, prtlevel = 1)

    if(wolfe_cond_res$fail == 1){break} # stop when reach minium in pk
    alpha_k = wolfe_cond_res$alpha
    xk = wolfe_cond_res$x      # update xk
    gk = wolfe_cond_res$grd    # update gk
    f_up = wolfe_cond_res$f
    fn_ls[[k]] = f_up
    x_list[[k]] = xk
    
    cat("\n", "At iterate", k, "\t", "f_up:", f_up)
    
    steps[k] = wolfe_cond_res$steps
    conds[k] = wolfe_cond_res$fail
    
    
    #orthog()
    uk = crossprod(Zk, gk)
    rho_z = gk - Zk %*% uk
    rhok <- sqrt(sum(rho_z^2))
    gk_1_accept = (rhok > accept_tol && rk < D)
    if (gk_1_accept){
      #orthog()
      rk = rk + 1
      zk = rho_z / rhok
      Zk = cbind(Zk, zk)
      Tk = cbind(Tk, uk)
      Tk = rbind(Tk, c(rep(0.0, rk - 1), rhok))
    } #else {}
    
    # expand()
    if(gk_1_accept){ #if rk+1 = rk + 1 #gk+1 is accept
      cat("\t", "gk+1 accept", "\t")
      Rk_prime = as.matrix( bdiag(Rk, sqrt(sigmak)))
      vk_prime = c(vk, 0.0)
      uk_prime = c(uk, rhok)
      qk_prime = c(qk, 0.0)
    } else {
      Rk_prime = Rk
      vk_prime = vk
      uk_prime = uk
      qk_prime = qk
    }
    
    #update()
    sk = alpha_k * qk_prime
    yk = uk_prime - vk_prime
    Lcs = Rk_prime %*% sk
    sum_yTs = sum(yk*sk)
    sum_Lcs2 = sum(Lcs^2)
    Hcs = crossprod(Rk_prime, Lcs)
    #w = (yk - sqrt(sum_yTs / sum_Lcs2) * (Rk_prime %*% Lcs))
    #z = Lcs / (sqrt(sum_yTs * sum_Lcs2))
    # check
    # Hk = crossprod(Rk_prime)
    # Hk_1 = Hk - (Hk%*%sk%*%t(sk)%*%Hk) /(t(sk)%*%Hk%*%sk)[1, 1] + 
    #   yk %*% t(yk) / sum(yk*sk)
    # cholHk_1 = chol(Hk_1)
    
    Rk_prime_prime = cholup(Rk_prime, yk / sqrt(sum_yTs), up = TRUE)
    Rk_prime_prime = cholup(Rk_prime_prime , Hcs / sqrt(sum_Lcs2), up = FALSE)
    # might find better function for the above step
    
    # reinitialize()
    sigmak = sum(yk^2) / sum_yTs    # setting follows Section 5.2
    Rk = Rk_prime_prime
    Rk[rk, rk] = sqrt(sigmak)
    vk = uk_prime
    
    # check convergence condition
    if(abs(f_up - f_old) < abs_tol){break}
    f_old = f_up
  }
  return(list(x_list = x_list, fn_ls = fn_ls, Rk = Rk, Zk = Zk, rank = rk, 
              steps = steps, 
         conds = conds, n = D, sigmak = sigmak, f = f_up))
}

RHR_monte_carlo <- function(lp_f,  # target log
                            fn,    # the minimize function
                            D,     # No. of pars
                            rk,    # rank of reduces hessian
                            Rk,    # chol from RHR
                            Zk,    # Zk from RHR
                            xk,    # xk from RHR
                            sigmak, # sigmak from RHR 
                            n_sam){
  
  v = matrix(rnorm(D * n_sam), nrow = D)
  u1 = crossprod(Zk, v)
  samples = Zk %*% backsolve(Rk, u1) + 1 / sqrt(sigmak) * (v - Zk %*% u1)
  samples = samples + xk[, 1]  # center at minimum found by RHR
  lp_list = apply(samples, 2, lp_f)
  fn_list = apply(samples, 2, fn)
  fn_min = fn(xk)
  E_lp = mean(lp_list); 
  E = mean(fn_list);
  next_flag = (min(fn_list) < fn_min)
  if(next_flag == FALSE){
    return(list(next_flag = next_flag, E_lp = E_lp, E = E))
  } else {
    return(list(next_flag = next_flag, E_lp = E_lp, E = E,
                x0_up = samples[, which.min(fn_list)]))
  }
} 


###--- test --- ###
# library('rstan')
# program <-
#   "
# data {
#   int D;
# } parameters {
#   vector[D] theta;
#   real<lower = 0.01> sigma;
# } model {
#   theta ~ normal(0, sigma);
#   sigma ~ normal(0, 1);
# }
# "
# D1 <- 40
# model <- stan_model(model_code = program)
# data = list(D = D1)
# 
# posterior <- to_posterior(model, data)
# D <- get_num_upars(posterior)
# init_bound = 2.0
# x0 <- runif(D, -init_bound, init_bound)
# fn <- function(theta) -log_prob(posterior, theta, adjust_transform = FALSE,
#                                 gradient = TRUE)[1]
# gr <- function(theta) -grad_log_prob(posterior, theta, adjust_transform = FALSE)
# lp_f <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE,
#                                  gradient = TRUE)[1]
# 
# sigma0 = 1.0
# max_iter = 300
# abs_tol = 1e-3
# accept_tol = 1e-18
# RHR_test <- RHR(fn, gr, x0, sigma0, max_iter, abs_tol, accept_tol)
# 
# tt <- optim(par = x0,
#             fn = fn,  # negate for maximization
#             gr = gr,
#             method = "L-BFGS-B",
#             control = list(maxit = 300, factr = 1e10, #ndeps = 1e-8 #,
#                            trace = 3, REPORT = 1))
# tt$value
# RHR_test$f
# 
# D = RHR_test$n
# rk = RHR_test$rank
# Rk = RHR_test$Rk
# Zk = RHR_test$Zk
# xk = RHR_test$x_list[[length(RHR_test$x_list)]]
# sigmak = RHR_test$sigmak
# n_sam = 100
# update <- RHR_monte_carlo(lp_f, fn, D, rk, Rk, Zk, xk, sigmak, n_sam)
# update


