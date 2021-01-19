library(loo)
library(Matrix)
printf <- function(msg, ...) cat(sprintf(msg, ...), "\n")

to_posterior <- function(model, data) {
  sampling(model, data = data, chains = 1, iter = 1, refresh = 0,
           algorithm = "Fixed_param")
}

# Returns optimization path for initialization x_init using objective
# function fn with gradient function gr, for a total of N iterations.
# Uses the L-BFGS-B algorithm of Nocedal.  Each row in the returned
# matrix is an iteration, and there are N + 1 rows because the
# initialization is included.
#
#
# @param init  initial parameter values
# @param fn          negative of log-density 
# @param gr          gradient function of negative log-density 
# @param lp_f        log-density
# @param init_bound  the boundwith of random initials for each dimension
# @param N1          maxium number of iterations for L-BFGS (default = 300)
# @param N_mode_max  maximum number of loop in pathfinder (default = 20)
# @param factr_tol   the option factr_tol in optim() (default = 1e9)
# @param lmm         the option lmm in optim() (default = 5)
# @return optimization path, matrix with N rows, each of which is a
#         point on the optimization path followed by objective
# E_lp = E_lp, E = E, lVol = lVol, log_MASS_c = log_MASS_c, 
# pareto_k = pareto_k, cond_num = cond_num, 
# y = y, step_count = step_count, fn_call = fn_call, 
# gr_call = gr_call, lgnorms = lgnorms, 
# sknorm_ls = sknorm_ls, thetak_ls = thetak_ls


opt_path <- function(init, fn, gr, 
                     init_bound = 2.0,
                     N1 = 300, N_mode_max = 20, 
                     N_sam = 100, factr_tol = 1e9, lmm = 5) {
  D <- length(init)
  #init <- runif(D, -init_bound, init_bound)
  # preallocation
  E_lp <- c()
  lVol <- c()
  log_MASS_c <- c()
  step_count = c()
  fn_call = c()     # No. calls to fn
  gr_call = c()     # No. calls to gr
  pareto_k  = c()
  lgnorms = c()     # log of gradient's norm
  sknorm_ls = c()    # distance of updates
  thetak_ls = c()     # curvature of updates
  
  
  y <- matrix(NA, nrow = 1, ncol = D + 1)
  y[1, 1:D] <- init
  y[1, D + 1] <- -fn(init)
  lgnorms <- c(lgnorms, log(sqrt(sum(gr(init)^2))))
  n0 = 1
  for (j in 1:N_mode_max){
    cat(j, "\t")
    if(j == 1){
      LBFGS_fail <- FALSE
      tryCatch(
        my_data <- capture.output(
          tt <- optim(par = y[n0, 1:D],
                      fn = fn,  # negate for maximization
                      gr = gr,
                      method = "L-BFGS-B",
                      control = list(maxit = N1, factr = factr_tol,
                                     pgtol = 0.0,
                                     #ndeps = 1e-8 #,
                                     trace = 6, REPORT = 1, lmm = lmm))
        ), 
         error = function(e) { LBFGS_fail <<- TRUE})
      if(LBFGS_fail){
        for(l in 1:40){
          print("\n reinitialize \n")
          y[n0, 1:D] <- runif(D, -init_bound, init_bound)
          y[n0, D + 1] <- -fn(y[n0, 1:D])
          lgnorms <- c(lgnorms[-length(lgnorms)], 
                       log(sqrt(sum(gr(y[n0, 1:D])^2))))
          LBFGS_fail <- FALSE
          tryCatch(my_data <- capture.output(
            tt <- optim(par = y[n0, 1:D],
                        fn = fn,  # negate for maximization
                        gr = gr,
                        method = "L-BFGS-B",
                        control = list(maxit = N1, factr = factr_tol, #ndeps = 1e-8 #,
                                       pgtol = 0.0,
                                       trace = 6, REPORT = 1, lmm = lmm))
          ), error = function(e) { LBFGS_fail <<- TRUE})
          if(!LBFGS_fail){break}
        }
      }
    } else { #j > 1
      LBFGS_fail <- FALSE
      tryCatch(my_data <- capture.output(
        tt <- optim(par = y[n0, 1:D],
                    fn = fn,  # negate for maximization
                    gr = gr,
                    method = "L-BFGS-B",
                    control = list(maxit = N1, factr = factr_tol, #ndeps = 1e-8 #,
                                   pgtol = 0.0, 
                                   trace = 6, REPORT = 1, lmm = lmm))
      ), error = function(e) { LBFGS_fail <<- TRUE})
    }
    
    if(LBFGS_fail){
      if(j == 1){
        E_lp[1] <- NA
        y = y
        return(list(E_lp = E_lp, lVol = lVol, log_MASS_c = log_MASS_c, 
                    pareto_k = pareto_k,  
                    y = y, step_count = step_count, fn_call = fn_call, 
                    gr_call = gr_call, lgnorms = lgnorms, 
                    sknorm_ls = sknorm_ls, thetak_ls = thetak_ls))
      }
      break} # fail in j > 1, break
    
    
    # recover the optimization trajectory and gradient
    L = length(my_data); L
    splited_output = unlist(lapply(my_data, f <- function(x){
      strsplit(as.character(x),split = " ")}))
    
    G_ind = which(splited_output == "G") + 2
    Iter = length(G_ind);
    if(tt$convergence == 0){Iter = Iter - 1}
    X = matrix(NA, nrow = Iter, ncol = D)
    G = matrix(NA, nrow = Iter, ncol = D)
    for(g in 1:Iter){
      X[g, ] = as.numeric(splited_output[(G_ind[g] - D - 2):(G_ind[g] - 3)])
      G[g, ] = as.numeric(splited_output[G_ind[g]:(G_ind[g] + D - 1)])
    }
    
    fn_ls <- apply(X, 1, fn)
    
    ### record geoinfo 
    Ykt = as.matrix(G[-1, ] - G[-nrow(G), ])
    Skt = as.matrix(X[-1, ] - X[-nrow(X), ])
    Dk = c()
    thetak = c()
    for(i in 1:(Iter - 1)){
      Dk[i] = sum(Ykt[i, ] * Skt[i, ])
      thetak[i] = sum(Ykt[i, ]^2) / Dk[i]   # curvature checking
    }
    
    sknorm_ls = c(sknorm_ls, sqrt(rowSums(Skt^2)))       
    thetak_ls = c(thetak_ls, thetak)   
    
    # find center of Gaussian
    lgnorm = log(sqrt(rowSums(G^2)))
    center_ind = which.min(lgnorm)   # center at points with smallest gradient norm
    lgnorms = c(lgnorms, lgnorm)
    
    
    # save results
    fn_ls <- apply(X, 1, fn)
    y = rbind(y, cbind(X, -fn_ls))
    step_count[j] <- tt$counts[1]
    fn_call[j] <- tt$counts[1]
    gr_call[j] <- tt$counts[2]
    
    
    if(center_ind < 3){
      print("on the edge of funnel")
      E_lp[j] = -fn_ls[center_ind];
      lVol[j] = -Inf
      log_MASS_c[j] = -Inf
      pareto_k[j] = -Inf
      
      return(list(E_lp = E_lp, lVol = lVol, log_MASS_c = log_MASS_c, 
                  pareto_k = pareto_k, 
                  y = y, step_count = step_count, fn_call = fn_call, 
                  gr_call = gr_call, lgnorms = lgnorms, 
                  sknorm_ls = sknorm_ls, thetak_ls = thetak_ls))
    }
    
    Ykt = as.matrix(G[2:center_ind, ] - G[1:(center_ind-1), ])
    Skt = as.matrix(X[2:center_ind, ] - X[1:(center_ind-1), ])
    Dk = c()
    thetak = c()
    m = center_ind - 1
    for(i in 1:m){
      Dk[i] = sum(Ykt[i, ] * Skt[i, ])
      thetak[i] = sum(Ykt[i, ]^2) / Dk[i]   # curvature checking
    }
    
    Neg_D_ind <- c()
    #curvature condition: abs(thetak) > 1/.Machine$double.eps
    if(any(Dk < 0 | abs(thetak) > 1e12)){   #any(Dk < 0 | abs(thetak) > 1e7)
      print("Negative or unstable Hessian")
      Neg_D_ind = which(Dk < 0 | abs(thetak) > 1e12)
      Ykt = as.matrix(Ykt[-Neg_D_ind, ])
      Skt = as.matrix(Skt[-Neg_D_ind, ])
      Dk = Dk[-Neg_D_ind]
      thetak = thetak[-Neg_D_ind]
      fn_ls_copy = fn_ls[1:center_ind][-(Neg_D_ind + 1)]
      m = length(fn_ls_copy) -
        which(fn_ls_copy < (fn_ls_copy[length(fn_ls_copy)] + 2*D ))[1] 
      m = min(Iter - 1, m)
      m = max(m, 2)
      x_center = as.matrix(X[((1:center_ind)[-(Neg_D_ind + 1)]), ])[length(fn_ls_copy), ]
      fn_center = fn_ls_copy[length(fn_ls_copy)]
    } else {
      m = center_ind - which(fn_ls[1:center_ind] < (fn_ls[center_ind] + 2*D))[1] 
      m = min(center_ind - 1, m)
      m = max(m, 2)
      x_center = X[center_ind, ]
      fn_center = fn_ls[center_ind]
    }
    
    
    # simulation samples of approximated posterior distribution and estimate E_lp
    ## build hessian ##
    lDk = length(Dk)
    #init_theta_ind = max(lDk - m, 1)
    #theta = sum(Ykt[init_theta_ind, ]^2) / 
    #  sum(Ykt[init_theta_ind, ]*Skt[init_theta_ind, ])
    Ykt = as.matrix(Ykt[(lDk - m + 1):lDk, ])
    Skt = as.matrix(Skt[(lDk - m + 1):lDk, ])
    Dk = Dk[(lDk - m + 1):lDk]
    thetak = thetak[(lDk - m + 1):lDk]
    
    ## take off sharp updates ##
    theta_sm <- c()
    small_lsk_ind <- c()
    if(m >= 10){ # if m >= 10 # m > lskupn
      lsk <- sqrt(rowSums(Skt^2))
      gm_lsk <- exp(mean(log(lsk))) / 4 #exp(mean(log(Dk))) #1/(mean(1/Dk))
      small_lsk_ind = which(lsk < gm_lsk)
      #small_lsk_ind = order(lsk)[(lskupn + 1):length(lsk)]
      if(length(small_lsk_ind) > 0){
        Ykt = as.matrix(Ykt[-small_lsk_ind, ])
        Skt = as.matrix(Skt[-small_lsk_ind, ])
        Dk = Dk[-small_lsk_ind]
        theta_sm <- thetak[small_lsk_ind]
        thetak = thetak[-small_lsk_ind]
        m = m - length(small_lsk_ind)
      }
    }
    
    # inv_theta_D = 1 / thetak[1] + Skt[1, ]^2 / Dk[1] - 
    #   Dk[1] * Ykt[1, ]^2 / sum(Ykt[1, ]^2)^2
    inv_theta_D = (thetak[1] + Ykt[1, ]^2 / Dk[1] -
                 thetak[1] * Skt[1, ]^2 / sum(Skt[1, ]^2))^{-1}
    for(d in 2:m){
      # inv_theta_D = Dk[d] * inv_theta_D / sum(Ykt[d, ]^2 * inv_theta_D) + 
      #   Skt[d, ]^2 / Dk[d] - 
      #   Dk[d] * (inv_theta_D *  Ykt[d, ])^2 / sum(inv_theta_D * Ykt[d, ]^2)^2
      inv_theta_D = 1 / ((sum(inv_theta_D * Ykt[d, ]^2) / Dk[d]) /
        inv_theta_D + Ykt[d, ]^2 / Dk[d] -
        (sum(inv_theta_D * Ykt[d, ]^2) / Dk[d]) * (Skt[d, ] / inv_theta_D)^2 /
        sum(Skt[d, ]^2 / inv_theta_D))
    }
    
    theta_D = 1 / inv_theta_D

    Rk = matrix(0.0, nrow = m, ncol = m)
    for(s in 1:m){
      for(i in 1:s){
        Rk[i, s] = sum(Skt[i, ] * Ykt[s, ])
      }
    }
    ninvRST = -backsolve(Rk, Skt)
    
    if( 2*m >= D){
      # directly calculate inverse Hessian and the cholesky decomposition
      Hk = diag(c(inv_theta_D), nrow = D) + 
        crossprod(Ykt %*% diag(inv_theta_D, nrow = D), ninvRST) + 
        crossprod(ninvRST, Ykt %*% diag(inv_theta_D, nrow = D))  + 
        crossprod(ninvRST, 
                  (diag(Dk) + 
                     tcrossprod(Ykt %*% diag(sqrt(inv_theta_D), nrow = D))) %*% 
                    ninvRST)
      cholHk = chol(Hk)
      logdetcholHk = determinant(cholHk)$modulus
      u = matrix(rnorm(D * N_sam), ncol = N_sam)
      u2 = crossprod(cholHk, u) + x_center
      
      # generate the samples uniformly from the region x^T H^{-1} x < qchisq(0.95, D)
      gamma = qchisq(0.95, D)
      s_V1 = sqrt(colSums(u^2))
      u_V1 = sapply(1:N_sam, function(i){ u[, i] / s_V1[i]})
      r_V1 = runif(N_sam)^{1/D}
      y_V = sapply(1:N_sam, function(i){ u_V1[, i] / r_V1[i]}) 
      z2 = sqrt(gamma)*crossprod(cholHk, y_V) + x_center
      
    } else {
       # use equation ?? to sample
      Wkbart = rbind(Ykt %*% diag(sqrt(inv_theta_D)), 
                     ninvRST %*% diag(sqrt(theta_D)))
      Mkbar = rbind(cbind(matrix(0.0, nrow = m, ncol = m), diag(m)),
                    cbind(diag(m), 
                          (diag(Dk) + 
                             tcrossprod(Ykt %*% diag(sqrt(inv_theta_D))))))
      qrW = qr(t(Wkbart))
      Qk = qr.Q(qrW)
      Rkbar = qr.R(qrW)
      Rktilde = chol(Rkbar %*% Mkbar %*% t(Rkbar) + diag(nrow(Rkbar)))
      logdetcholHk = sum(log(diag(Rktilde))) - 0.5 * sum(log(theta_D))
      u = matrix(rnorm(D * N_sam), ncol = N_sam)
      u1 = crossprod(Qk, u)
      u2 = diag(sqrt(inv_theta_D)) %*% 
        (Qk %*% crossprod(Rktilde, u1) + (u - Qk %*% u1)) + x_center
      
      # generate the samples uniformly from the region x^T H^{-1} x < qchisq(0.95, D)
      gamma = qchisq(0.95, D)
      s_V1 = sqrt(colSums(u^2))
      u_V1 = sapply(1:N_sam, function(i){ u[, i] / s_V1[i]})
      r_V1 = runif(N_sam)^{1/D}
      y_V = sapply(1:N_sam, function(i){ u_V1[, i] / r_V1[i]}) 
      z1 = crossprod(Qk, y_V)
      z2 = diag(sqrt(gamma)/sqrt(theta_D)) %*% 
        (Qk %*% crossprod(Rktilde, z1) + (y_V - Qk %*% z1)) + x_center
      
    }
    
    fn_draws <-  c()
    ind_draws <- c()
    lp_approx_draws <- c()
    lp_f_z2_draws <- c() # records the log_density of samples uniformly from
                         # the 95% CI of approximating Gaussian
    for(l in 1:ncol(u)){
      
      # skip bad samples
      skip_flag = FALSE
      tryCatch(f_test <- fn(u2[, l]),
               error = function(e) { skip_flag <<- TRUE})
      if(skip_flag){next} 
      else {
        
        ind_draws <- c(ind_draws, l)
        fn_draws <- c(fn_draws, f_test)
        lp_approx_draws <- c(lp_approx_draws, - logdetcholHk - 
                               0.5 * sum(u[, l]^2))
        lp_f_z2_draws <- c(lp_f_z2_draws, -fn(z2[, l]))
      }
    }
    fn_call[j] = fn_call[j] + N_sam
    ### psis test ###
    if(is.infinite(max(fn_draws))){
      keep_ind = which(is.finite(fn_draws))
      fn_draws = fn_draws[keep_ind]
      lp_approx_draws = lp_approx_draws[keep_ind]
    } 
    lp_ratios = - fn_draws - lp_approx_draws 
    psis_test = psis(lp_ratios, r_eff = 1)
    pareto_k[j] = pareto_k_values(psis_test)

    E_lp[j] = sum(c(weights(psis_test, log = FALSE)) * -fn_draws) #mean(lp_f_draws)

    lVol[j] = logdetcholHk #determinant(cholHk)$modulus #-0.5*determinant(Bk)$modulus
    # log_MASS_c[j] = lVol[j] + log(sum(c(unif_weight) * 
    #                                     (exp(-fn_draws + mean(fn_draws))))) - 
    #   mean(fn_draws)
    log_MASS_c[j] = lVol[j] + 
      log(mean(exp(lp_f_z2_draws[is.finite(lp_f_z2_draws)] - 
                     mean(lp_f_z2_draws[is.finite(lp_f_z2_draws)])))) + 
      mean(lp_f_z2_draws[is.finite(lp_f_z2_draws)])
    
    if(is.na(log_MASS_c[j]) | is.infinite(log_MASS_c[j])){log_MASS_c[j] = -Inf}
    
    if(min(fn_draws) < fn_center && (j < N_mode_max)){  #tt$value
      n0 <- nrow(y) + 1
      init <- u2[, which.min(fn_draws)]
      y = rbind(y, c(init, -fn(init)))
      lgnorms = c(lgnorms, log(sum(gr(init)^2)))
    } else {
      break
    }
  }
  #E_lp
  return(list(E_lp = E_lp, lVol = lVol, log_MASS_c = log_MASS_c, 
              pareto_k = pareto_k,  
              y = y, step_count = step_count, fn_call = fn_call, 
              gr_call = gr_call, lgnorms = lgnorms, 
              sknorm_ls = sknorm_ls, thetak_ls = thetak_ls))
}


# Returns optimization path for specified Stan model and data for
# the specified number of iterations using the specified bound on
# uniform initialization on the unconstrained scale.  See opt_path()
# for a description of output format and algorithm.
#
# @param model Stan model (compiled using rstan::stan_model())
# @param data data list for call to rstan::sampling for specified model
# @param N total number of iterations (default = 25)
# @param init_bound upper bound on uniform initialization; negation is lower
# @return optimization path (see opt_path documentation)
opt_path_stan <- function(model, data, N1, N_mode_max, N_sam, 
                          factr_tol = 1e9, lmm = 5,
                          init_bound = 2) {
  # require chains > 0, iter > 0, so use Fixed_param to avoid work
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                 gradient = TRUE)[1] 
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  out <- opt_path(init, fn = fn, gr = gr, N1 = N1,
                  N_mode_max = N_mode_max, N_sam = N_sam, factr_tol = factr_tol,
                  lmm = lmm)
  return(out)
}

#opt_path_stan(model, data, N1, N_rep, init_bound = 2)

opt_path_stan_parallel <- function(seed_list, mc.cores, model, data, 
                                   N1, N_mode_max, N_sam, init_bound, 
                                   factr_tol, lmm){
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)

  MC = length(seed_list)
  init = c()
  for(i in 1:MC){
    set.seed(seed_list[i])
    init[[i]] <- runif(D, -init_bound, init_bound)
  }
  out <- mclapply(init, opt_path, fn = fn, gr = gr, N1 = N1,
                  N_mode_max = N_mode_max, N_sam = N_sam, factr_tol = factr_tol,
                  lmm = lmm, mc.cores = mc.cores)
}

opt_path_stan_init_parallel <- function(init_ls, mc.cores, model, data, 
                                        N1, N_mode_max, N_sam, init_bound, 
                                        factr_tol, lmm){
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  out <- mclapply(init_ls, opt_path, fn = fn, gr = gr, N1 = N1,
                  N_mode_max = N_mode_max, N_sam = N_sam, factr_tol = factr_tol,
                  lmm = lmm, mc.cores = mc.cores)
}


# Return optimization path with last column (objective function value)
# removed.
#
# @param path optimization path with last column for objective
#        function value
# @return
params_only <- function(path) {
  N <- dim(path)[1]
  D <- dim(path)[2]
  path[1:N, 1:(D - 1)]
}

find_indx <- function(param_path) {
  mode_ind = which.max(param_path$log_MASS_c)
  E_target = param_path$E_lp[mode_ind]
  mark_index = which.min(abs(param_path$y[, (ncol(param_path$y))] - 
                                param_path$E_lp[mode_ind]))
  return(mark_index)
}

fit_info <- function(param_path) {
  # extract fitting information #
  
  mode_ind = which.max(param_path$log_MASS_c)
  E_target = param_path$E_lp[mode_ind]
  mark_index = which.min(abs(param_path$y[, (ncol(param_path$y))] - E_target))
  log_MASS_c = param_path$log_MASS_c[mode_ind]
  pareto_k = param_path$pareto_k[mode_ind]
  
  return(c(log_MASS_c, pareto_k, E_target))
}



# library('rstan')
# program <-
#   "
# data {
#   int D;
# } parameters {
#   vector[D] theta;
#   real<lower = 0> sigma;
# } model {
#   theta ~ normal(0, sigma);
#   sigma ~ normal(0, 1);
# }
# "
# D <- 100
# model <- stan_model(model_code = program)
# data = list(D = D)
# opath <- opt_path_stan(model, data = data, N = 75, init_bound = 5)
# find_typical(model, data, opath)

