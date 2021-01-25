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
                     N_sam = 1, N_mass = 5, factr_tol = 1e9, lmm = 5) {
  D <- length(init)
  #init <- runif(D, -init_bound, init_bound)
  # preallocation
  E_lp <- c()
  lVol <- c()
  log_MASS_c <- c()
  ELBO <- c()       # record ELBO
  sample_pkg_save <- list() # save the sample info for each mode
  ELBO_save <- list() # save the ELBO related info for each mode
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
        y = y
        return(list(log_MASS_c = log_MASS_c, sample_pkg_save = sample_pkg_save, 
                    ELBO_save = ELBO_save, y = y, fn_call = fn_call, 
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
    
    # save results
    lgnorm = log(sqrt(rowSums(G^2)))
    lgnorms = c(lgnorms, lgnorm)
    fn_ls <- apply(X, 1, fn)
    y = rbind(y, cbind(X, -fn_ls))
    fn_call[j] <- tt$counts[1]
    gr_call[j] <- tt$counts[2]
    fn_draws <-  c()
    sample_draws <- rep(0.0, D)
    
    # estimate ELBO and save results for the maximum ELBO
    ELBO_ls <-  c()
    ELBO_max <- -Inf
    sample_pkg_pick <- list()
    t <- proc.time()
    for (iter in 2:(Iter-1)){
      cat(iter, "\t")
      ill_distr = FALSE
      tryCatch(sample_pkg <- 
                 Get_H(X[iter+1, ], Ykt[1:iter, ], Skt[1:iter, ],
                       Dk[1:iter], thetak[1:iter], sknorm_ls[1:iter], 
                       thetak_ls[1:iter], fn_ls[1:(iter + 1)]),
               error = function(e) { ill_distr <<- TRUE})
      if(ill_distr){ next }
      if(is.na(sample_pkg[1])){ next }
      ELBO_fit <- est_ELBO(sample_pkg, N_sam, D, fn)
      fn_call[j] <- fn_call[j] + ELBO_fit$fn_calls_ELBO
      fn_draws <- c(fn_draws, ELBO_fit$fn_draws)
      sample_draws <- cbind(sample_draws, ELBO_fit$repeat_draws)
      ELBO_ls <- c(ELBO_ls, ELBO_fit$ELBO)
      if(ELBO_fit$ELBO > ELBO_max){
        ELBO_max <- ELBO_fit$ELBO
        ELBO_fit_pick <- ELBO_fit
        sample_pkg_pick <- sample_pkg
      }
    }
    proc.time() -t
    sample_pkg_save[[j]] <- sample_pkg_pick
    ELBO_save[[j]] <- ELBO_fit_pick
    plot((1:length(ELBO_ls))+2, ELBO_ls, type = "l", ylim = c(-10, 20), 
         xlab = "iter", ylab = "ELBO")
    # print(-ELBO_fit_pick$fn_draws)  # 2661
    # lp_mean[1]
    
    # estimate log probability mass
    log_MASS_fit <- est_log_mass(sample_pkg_pick, N_mass, D, fn)
    log_MASS_c[j] <- log_MASS_fit$log_MASS_c
    fn_call[j] <- fn_call[j] + log_MASS_fit$fn_calls_mass
    
    # check for next update
    if(min(fn_draws) < tt$value && (j < N_mode_max)){  #tt$value
      n0 <- nrow(y) + 1
      init <- sample_draws[, which.min(fn_draws) + 1]
      y = rbind(y, c(init, -fn(init)))
      lgnorms = c(lgnorms, log(sum(gr(init)^2)))
    } else {
      break
    }
  }
  #E_lp
  return(list(log_MASS_c = log_MASS_c, sample_pkg_save = sample_pkg_save, 
              ELBO_save = ELBO_save, y = y, fn_call = fn_call, 
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
opt_path_stan <- function(model, data, N1, N_mode_max, N_sam, N_mass,
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
                                   N1, N_mode_max, N_sam, N_mass, init_bound, 
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
                  N_mode_max = N_mode_max, N_sam = N_sam, N_mass = N_mass, 
                  factr_tol = factr_tol,
                  lmm = lmm, mc.cores = mc.cores)
}

opt_path_stan_init_parallel <- function(init_ls, mc.cores, model, data, 
                                        N1, N_mode_max, N_sam, N_mass, 
                                        init_bound, factr_tol, lmm){
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  out <- mclapply(init_ls, opt_path, fn = fn, gr = gr, N1 = N1,
                  N_mode_max = N_mode_max, N_sam = N_sam, N_mass = N_mass,
                  factr_tol = factr_tol, lmm = lmm, mc.cores = mc.cores)
}


Get_H <- function(x_center, Ykt, Skt, Dk, thetak, sknorm_ls, thetak_ls, fn_ls){
  
  D = length(x_center)
  fn_pt = fn_ls[length(fn_ls)]
  fn_ls_copy = fn_ls
  
  Neg_D_ind <- c()
  #curvature condition: abs(thetak) > 1/.Machine$double.eps
  if(any(Dk < 0 | abs(thetak) > 1e12)){   #any(Dk < 0 | abs(thetak) > 1e7)
    print("Negative or unstable Hessian")
    Neg_D_ind = which(Dk < 0 | abs(thetak) > 1e12)
    Ykt = as.matrix(Ykt[-Neg_D_ind, ])
    Skt = as.matrix(Skt[-Neg_D_ind, ])
    Dk = Dk[-Neg_D_ind]
    thetak = thetak[-Neg_D_ind]
    fn_ls_copy = fn_ls[-(Neg_D_ind + 1)]
  }
  m = length(Dk)
  
  ######
  # take off log density too far from the center
  # m = length(fn_ls_copy) -
  #   which(fn_ls_copy < (fn_pt + 2*D ))[1]
  # 
  # lDk = length(Dk)
  # Ykt = as.matrix(Ykt[(lDk - m + 1):lDk, ])
  # Skt = as.matrix(Skt[(lDk - m + 1):lDk, ])
  # Dk = Dk[(lDk - m + 1):lDk]
  # thetak = thetak[(lDk - m + 1):lDk]
  # 
  # ## take off sharp updates ##
  # theta_sm <- c()
  # small_lsk_ind <- c()
  # if(m >= 10){ # if m >= 10 # m > lskupn
  #   lsk <- sqrt(rowSums(Skt^2))
  #   gm_lsk <- exp(mean(log(lsk))) / 5 #exp(mean(log(Dk))) #1/(mean(1/Dk))
  #   small_lsk_ind = which(lsk < gm_lsk)
  #   #small_lsk_ind = order(lsk)[(lskupn + 1):length(lsk)]
  #   if(length(small_lsk_ind) > 0){
  #     Ykt = as.matrix(Ykt[-small_lsk_ind, ])
  #     Skt = as.matrix(Skt[-small_lsk_ind, ])
  #     Dk = Dk[-small_lsk_ind]
  #     theta_sm <- thetak[small_lsk_ind]
  #     thetak = thetak[-small_lsk_ind]
  #     m = m - length(small_lsk_ind)
  #   }
  # }
  ######
  
  if(m < 2){ # cannot approximate Hessian
    return(NA)
  }
  
  # simulation samples of approximated posterior distribution and estimate E_lp
  inv_theta_D = (thetak[1] + Ykt[1, ]^2 / Dk[1] -
                   thetak[1] * Skt[1, ]^2 / sum(Skt[1, ]^2))^{-1}
  for(d in 2:m){
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
    
    sample_pkg <- list(label = "full", cholHk = cholHk, 
                       logdetcholHk = logdetcholHk,
                       x_center = x_center)
    
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
    
    sample_pkg <- list(label = "sparse", theta_D = theta_D, 
                       Qk = Qk, Rktilde = Rktilde,  
                       logdetcholHk = logdetcholHk,
                       Mkbar = Mkbar,
                       Wkbart = Wkbart, 
                       x_center = x_center)
    
  }
  
  return(sample_pkg)
}

est_ELBO <- function(sample_pkg, N_sam, D, fn){
  
  fn_draws <-  rep(Inf, N_sam)
  lp_approx_draws <- rep(0.0, N_sam)
  repeat_draws <- matrix(0, nrow = D, ncol = N_sam)
  draw_ind = 1
  fn_calls_ELBO = 0
  f_test_ELBO = -Inf
  
  for(l in 1:(2*N_sam)){
    if(sample_pkg$label == "full"){
      u = rnorm(D)
      u2 = crossprod(sample_pkg$cholHk, u) + sample_pkg$x_center
    }else{
      u = rnorm(D)
      u1 = crossprod(sample_pkg$Qk, u)
      u2 = diag(sqrt(1 / sample_pkg$theta_D)) %*%
        (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, u1) + 
           (u - sample_pkg$Qk %*% u1)) + sample_pkg$x_center
    }
    # skip bad samples
    skip_flag = FALSE
    tryCatch(f_test_ELBO <- fn(u2),
             error = function(e) { skip_flag <<- TRUE})
    if(is.nan(f_test_ELBO)){skip_flag <<- TRUE}
    if(skip_flag){
      next
      } else {
      fn_draws[draw_ind] <- f_test_ELBO
      lp_approx_draws[draw_ind] <- - sample_pkg$logdetcholHk - 0.5 * sum(u^2)
      repeat_draws[, draw_ind] <- u2
      draw_ind = draw_ind + 1
    }
    fn_calls_ELBO = fn_calls_ELBO + 1
    if(draw_ind == N_sam + 1){break}
  }
  
  ### ELBO estimate ###
  #ELBO <- -mean(fn_draws) + sample_pkg$logdetcholHk + D/2
  ELBO <- -mean(fn_draws) - mean(lp_approx_draws)
  
  if(is.nan(ELBO)){ELBO <- -Inf}
  if(any(is.nan(fn_draws))){
    fn_draws[is.nan(fn_draws)] <- Inf
  }
  
  return(list(ELBO = ELBO, repeat_draws = repeat_draws, fn_draws = fn_draws, 
              lp_approx_draws = lp_approx_draws, fn_calls_ELBO = fn_calls_ELBO))
}

est_log_mass <- function(sample_pkg, N_sam, D, fn){
  gamma = qchisq(0.95, D)
  draw_ind = 1
  fn_calls_mass = 0
  lp_f_draws_mass <-  rep(-Inf, N_sam)
  
  for(l in 1:(2*N_sam)){
    
    # generate the samples uniformly from the region x^T H^{-1} x < qchisq(0.95, D)
    u = rnorm(D)
    s_V1 = sqrt(sum(u^2))
    u_V1 = u / s_V1
    r_V1 = runif(1)^{1/D}
    y_V = u_V1 / r_V1
    
    if(sample_pkg$label == "full"){
      z2 = sqrt(gamma) * crossprod(sample_pkg$cholHk, y_V) + sample_pkg$x_center
    }else{
      z1 = crossprod(sample_pkg$Qk, y_V)
      z2 = diag(sqrt(gamma)/sqrt(sample_pkg$theta_D)) %*%
        (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, z1) + 
           (y_V - sample_pkg$Qk %*% z1)) + sample_pkg$x_center
    }
    
    # skip bad samples
    skip_flag = FALSE
    tryCatch(f_test <- fn(z2),
             error = function(e) { skip_flag <<- TRUE})
    if(skip_flag){next} 
    else {
      lp_f_draws_mass[draw_ind] <- -f_test
      draw_ind = draw_ind + 1
    }
    fn_calls_mass = fn_calls_mass + 1
    if(draw_ind == N_sam + 1){break}
  }
  
  log_MASS_c = sample_pkg$logdetcholHk +
    log(mean(exp(lp_f_draws_mass[is.finite(lp_f_draws_mass)] -
                   mean(lp_f_draws_mass[is.finite(lp_f_draws_mass)])))) +
    mean(lp_f_draws_mass[is.finite(lp_f_draws_mass)])

  if(is.na(log_MASS_c) | is.infinite(log_MASS_c)){log_MASS_c = -Inf}
  
  return(list(log_MASS_c = log_MASS_c, fn_calls_mass = fn_calls_mass))
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
  log_MASS_c = param_path$log_MASS_c[mode_ind]

  return(c(log_MASS_c))
}

filter_mode <- function(param_path){
  fit_info_sum <- sapply(param_path, fit_info)
  pick_mode <- 1:length(fit_info_sum)
  if(any(is.finite(fit_info_sum)) & (sum(is.finite(fit_info_sum)) > 2)){
    finit_ind = which(is.finite(fit_info_sum))
    shapW_test_result <- shapiro.test(fit_info_sum[finit_ind])
    if(shapW_test_result$p.value < 0.01){
      kmean_results <- kmeans(fit_info_sum[finit_ind], centers = 2)
      k = which.max(kmean_results$centers)
      
      pick_mode = finit_ind[kmean_results$cluster == k]
      plot(finit_ind, fit_info_sum[finit_ind], col = kmean_results$cluster,
           type = "p", lwd = 2, ylab = "log scaled prob mass")
    }else{
      pick_mode = finit_ind
      plot(finit_ind, fit_info_sum[finit_ind],
           type = "p", lwd = 2, ylab = "log scaled prob mass")
    }
  }else{
    pick_mode = NA
    return(pick_mode)
  }
  
  for(lop in 1:(length(fit_info_sum) - 2)){
    if(length(pick_mode > 2)){
      shapW_test_result <- shapiro.test(fit_info_sum[pick_mode])
      if(shapW_test_result$p.value < 0.01){
        kmean_results <- kmeans(fit_info_sum[pick_mode], centers = 2)
        k = which.max(kmean_results$centers)
        plot(pick_mode, fit_info_sum[pick_mode], col = kmean_results$cluster,
             type = "p", lwd = 2, ylab = "log scaled prob mass")
        pick_mode = pick_mode[kmean_results$cluster == k]
      }else{break}
    }else{break}
  }
  return(pick_mode)
}

extract_samples <- function(param_path){
  mode_ind = which.max(param_path$log_MASS_c)
  sample <- param_path$ELBO_save[[mode_ind]]$repeat_draws
  return(sample)
}

extract_lps <- function(param_path){
  mode_ind = which.max(param_path$log_MASS_c)
  lps <- -param_path$ELBO_save[[mode_ind]]$fn_draws
  return(lps)
}


