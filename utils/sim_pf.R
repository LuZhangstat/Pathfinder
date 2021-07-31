library(loo)
library(Matrix)
library(matrixStats)
library(optimx)
printf <- function(msg, ...) cat(sprintf(msg, ...), "\n")

to_posterior <- function(model, data) {
  sampling(model, data = data, chains = 1, iter = 1, refresh = 0,
           algorithm = "Fixed_param")
}


opt_path <- function(init, fn, gr, 
                     init_bound = 2.0,
                     N1 = 1000,
                     N_sam_DIV = 5, N_sam = 100,
                     factr_tol = 1e2, lmm = 6, seed = 1234,
                     eval_lp_draws = TRUE) {
  
  #' Run one-path Pathfinder for initialization init 
  #' 
  #' @param init        initial parameter values
  #' @param fn          negative of log-density 
  #' @param gr          gradient function of negative log-density 
  #' @param init_bound  the boundwith of random initials for each dimension (default = 2.0)
  #' @param N1          maxium number of iterations for L-BFGS (default = 1000)
  #' @param N_sam_DIV   number of samples used in estimating ELBO (default = 5)
  #' @param N_sam       number of samples from the approximating Normal returned
  #'                    will not do extra samples when < N_sam_DIV  (default = 5)
  #' @param factr_tol   the option factr in optim() (default = 1e2)
  #' @param lmm         the option lmm in optim() (default = 6)
  #' @param seed        random seed for one-path Pathfinder
  #' @param eval_lp_draws logical; if TRUE (default) evaluate the log-densities of approximating draws
  #' @return 
  
  
  set.seed(seed)
  D <- length(init)
  # preallocation
  sample_pkg_save <- NA # save the sample info for each mode
  DIV_save <- -Inf # save the divergence related info for each mode
  fn_call = 0     # No. calls to fn
  gr_call = 0     # No. calls to gr
  
  ## retrieve the optimization trajectory from optim() ##
  y <- matrix(NA, nrow = 1, ncol = D + 1) # record the optimization trajectory and log-density
  LBFGS_fail <- TRUE
  N_try <- 0
  while(LBFGS_fail){   # if L-BFGS cannot run with the initial, find another inital
    LBFGS_fail <- FALSE
    if(N_try == 0){
      N_try = N_try + 1
      y[1, 1:D] <- init
    }else{
      print("\n reinitialize \n")
      N_try = N_try + 1
      y[1, 1:D] <- runif(D, -init_bound, init_bound)
    }
    tryCatch(y[1, D + 1] <- -fn(y[1, 1:D]), 
             error = function(e) { LBFGS_fail <<- TRUE})
    if(LBFGS_fail | is.infinite(y[1, D + 1])){
      LBFGS_fail <- TRUE
      fn_call = fn_call + 1  # record the evaluation of log-density of ill initials
      next}
    g1 <- gr(y[1, 1:D]) # record the gradient of initials
    if(any(is.na(g1))){
      gr_call = gr_call + 1  # record the evaluation of gradient of ill initials
      LBFGS_fail <- TRUE
      next
    }
    tryCatch(
      my_data <- capture.output(
        tt <- optimx(par = y[1, 1:D],
                     fn = fn,  # negate for maximization
                     gr = gr,
                     method = "L-BFGS-B",
                     control = list(maxit = N1, 
                                    pgtol = 0.0, 
                                    factr = factr_tol,
                                    trace = 6, REPORT = 1, lmm = lmm)),
        type = "output"), 
      error = function(e) { LBFGS_fail <<- TRUE})
    if(LBFGS_fail){ # fail the section of code that does the checking
      fn_call = fn_call + 1
      gr_call = gr_call + 1
      next}
    if(tt$convcode == 9999){
      LBFGS_fail <- TRUE
      splited_output = unlist(lapply(my_data, f <- function(x){
        strsplit(as.character(x),split = " ")}))
      
      line_search_ind = which(splited_output == "SEARCH") + 1 # check line seach times the 
      
      eval_count <- sum(as.numeric(splited_output[line_search_ind])) + 
        length(line_search_ind) + 1
      fn_call = fn_call + eval_count  
      gr_call = gr_call + eval_count
    }
  }
  
  # recover the optimization trajectory X and gradient G.
  L = length(my_data); L
  splited_output = unlist(lapply(my_data, f <- function(x){
    strsplit(as.character(x),split = " ")}))
  
  G_ind = which(splited_output == "G") + 2
  Iter = length(G_ind);
  if(tt$convcode == 0){Iter = Iter - 1}
  X = matrix(NA, nrow = Iter + 1, ncol = D)
  G = matrix(NA, nrow = Iter + 1, ncol = D)
  X[1, ] = y[1, 1:D]
  G[1, ] = g1  # add this since we cannot retrieve the gradient of the initial 
  for(g in 1:Iter){
    X[g + 1, ] = as.numeric(splited_output[(G_ind[g] - D - 2):(G_ind[g] - 3)])
    G[g + 1, ] = as.numeric(splited_output[G_ind[g]:(G_ind[g] + D - 1)])
  }
  
  ### record geometry info 
  Ykt = as.matrix(G[-1, ] - G[-nrow(G), ])
  Skt = as.matrix(X[-1, ] - X[-nrow(X), ])
  
  # save results
  fn_ls <- apply(X, 1, fn)         # record fn of the optimization trajectory
  y = rbind(y, cbind(X, -fn_ls))   # update y
  fn_call <- fn_call + tt$fevals          # record the calls of fn in L-BFGS
  gr_call <- gr_call + tt$gevals          # record the calls of gr in L-BFGS
  
  # estimate DIV for all approximating Gaussians and save results of the one with maximum DIV
  DIV_ls <-  c()
  DIV_fit_pick <- NA
  DIV_max <- -Inf
  sample_pkg_pick <- list()
  E <- rep(1, D)    # initial diag inv Hessian
  Ykt_h <- NULL; Skt_h <- NULL # initialize matrics for storing history of updates
  t <- proc.time()
  for (iter in 1:Iter){
    #cat(iter, "\t")
    inc_flag <- check_cond(Ykt[iter, ], Skt[iter, ])
    if(inc_flag){
      E <- Form_init_Diag(E, Ykt[iter, ], Skt[iter, ]) # initial estimate of diagonal inverse Hessian
      Ykt_h <- updateYS(Ykt_h = Ykt_h, Ykt = Ykt[iter, ], lmm) # update Y and S matrix
      Skt_h <- updateYS(Ykt_h = Skt_h, Ykt = Skt[iter, ], lmm)
    }
    ill_distr = FALSE
    tryCatch(
      # generate matrics for forming approximted inverse Hessian
      sample_pkg <- Form_N_apx_taylor(X[iter + 1, ], G[iter + 1, ], Ykt_h, Skt_h, E, lmm),
      #sample_pkg <- Form_N_apx(X[iter + 1, ], Ykt_h, Skt_h, E, lmm),
      error = function(e) { ill_distr <<- TRUE})
    if(ill_distr){ next }
    if(is.na(sample_pkg[1])){ next }
    DIV_label = "ELBO"; 
    DIV_fit <- est_DIV(sample_pkg, N_sam_DIV, fn, label = DIV_label) #lCDIV #lADIV  #lIKL #ELBO
    fn_call <- fn_call + DIV_fit$fn_calls_DIV
    DIV_ls <- c(DIV_ls, DIV_fit$DIV)
    if(DIV_fit$DIV > DIV_max){
      DIV_max <- DIV_fit$DIV
      DIV_fit_pick <- DIV_fit
      sample_pkg_pick <- sample_pkg
    }
  } 
  proc.time() -t 
  #plot(1:length(DIV_ls), DIV_ls, ylim = c(-15000, 2000))
  sample_pkg_save <- sample_pkg_pick
  DIV_save <- DIV_fit_pick
  
  if(is.null(DIV_ls) | (max(DIV_ls) == -Inf)){ # if no normal approximation generated, return the set of parameters found by optim
    return(list(sample_pkg_save = sample_pkg_save, # matrics for estiamted inv Hessian 
                DIV_save = DIV_save,               # approximating draws and DIV
                y = y,                             # optimization path
                fn_call = fn_call,                 # no. calls of log density
                gr_call = gr_call,                 # no. calls of gradient
                status = "mode"                    # status of the one-path pathfinder
    ))
  }
  
  ## Generate upto N_sam samples from the picked approximating Normal ##
  if(N_sam > length(DIV_fit_pick$fn_draws)){
    draws_N_apx <- Sam_N_apx(sample_pkg_save, 
                             N_sam - length(DIV_fit_pick$fn_draws))
    
    ## update the samples in DIV_save ##
    DIV_save$repeat_draws <- cbind(DIV_save$repeat_draws, draws_N_apx$samples)
    DIV_save$lp_approx_draws <- c(DIV_save$lp_approx_draws,
                                  draws_N_apx$lp_apx_draws)
    
    if(eval_lp_draws == TRUE){
      # if evaluate the log-density of approximating draws
      fn_draws_apx <- apply(draws_N_apx$samples, 2, 
                            f <- function(x){
                              ill_fn = FALSE
                              tryCatch(
                                nld <- fn(x),
                                error = function(e) { ill_fn <<- TRUE})
                              ifelse(ill_fn, Inf, nld)
                            })
      DIV_save$fn_draws <- c(DIV_save$fn_draws, fn_draws_apx)
      fn_call = fn_call +  N_sam - length(DIV_fit_pick$fn_draws)
    }
  }
  
  return(list(sample_pkg_save = sample_pkg_save, 
              DIV_save = DIV_save, 
              y = y, 
              fn_call = fn_call, 
              gr_call = gr_call,
              status = "samples"))
}


opt_path_stan <- function(model, data, init_bound = 2, N1 = 1000,
                          N_sam_DIV = 5, N_sam = 100,
                          factr_tol = 1e2, lmm = 6, seed = 1234,
                          eval_lp_draws = TRUE){
  
  #' Run one-path Pathfinder for specified Stan model and data for
  #' the specified number of iterations using the specified bound on
  #' uniform initialization on the unconstrained scale.  See opt_path()
  #' for a description of output format and algorithm.
  #'
  #' @param model          Stan model (compiled using rstan::stan_model())
  #' @param data           data list for call to rstan::sampling for specified model
  #' @param init_bound     the boundwith of random initials for each dimension (default = 2.0)
  #' @param N1             maxium number of iterations for L-BFGS (default = 1000)
  #' @param N_sam_DIV      number of samples used in estimating ELBO (default = 5)
  #' @param N_sam          number of samples from the approximating Normal returned
  #'                       will not do extra samples when < N_sam_DIV  (default = 5)
  #' @param factr_tol      the option factr in optim() (default = 1e2)
  #' @param lmm            the option lmm in optim() (default = 6)
  #' @param seed           random seed for generating initials
  #' @param eval_lp_draws  logical; if TRUE (default) evaluate the log-densities of approximating draws
  #' @return 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  set.seed(seed)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1] 
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  out <- opt_path(init, fn = fn, gr = gr, N1 = N1, N_sam_DIV = N_sam_DIV,
                  N_sam = N_sam, factr_tol = factr_tol,
                  lmm = lmm, seed = seed, eval_lp_draws = eval_lp_draws)
  #seed = seed + 1
  return(out)
}


opt_path_stan_parallel <- function(seed_init, seed_list, mc.cores, model, data, 
                                   init_bound = 2.0, N1 = 1000, 
                                   N_sam_DIV = 5, N_sam = 100, 
                                   factr_tol = 1e2, lmm = 6,
                                   eval_lp_draws = TRUE){
  
  #' Run one-path Pathfinder with random initials in parallel
  #'
  #' @param seed_init      array of random seeds for generating random initials
  #' @param seed_list      array of random seeds for running one-path Pathfinders
  #' @param mc.cores       The maximum number of one-path Pathfinder to run in parallel 
  #' @param model          Stan model (compiled using rstan::stan_model())
  #' @param data           data list for call to rstan::sampling for specified model
  #' @param init_bound     the boundwith of random initials for each dimension (default = 2.0)
  #' @param N1             maxium number of iterations for L-BFGS (default = 1000)
  #' @param N_sam_DIV      number of samples used in estimating ELBO (default = 5)
  #' @param N_sam          number of samples from the approximating Normal returned
  #'                       will not do extra samples when < N_sam_DIV  (default = 5)
  #' @param factr_tol      the option factr in optim() (default = 1e2)
  #' @param lmm            the option lmm in optim() (default = 6)
  #' @param eval_lp_draws  logical; if TRUE (default) evaluate the log-densities of approximating draws
  #' @return 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  MC = length(seed_init)
  init = c()
  list_ind = c()
  for(i in 1:MC){
    set.seed(seed_init[i])
    init[[i]] <- runif(D, -init_bound, init_bound)
    list_ind[[i]] <- i
  }
  
  out <- mclapply(list_ind, f <- function(x){
    opt_path(init = init[[x]] ,fn = fn, gr = gr, N1 = N1, N_sam_DIV = N_sam_DIV,
             N_sam = N_sam,  factr_tol = factr_tol,
             lmm = lmm, seed = seed_list[x], eval_lp_draws = eval_lp_draws)
  }, mc.cores = mc.cores)
}

opt_path_stan_init_parallel <- function(init_ls, mc.cores, model, data, 
                                        init_bound = 2.0, N1 = 1000, 
                                        N_sam_DIV = 5, N_sam = 100, 
                                        factr_tol = 1e2, lmm = 6,
                                        seed_list, eval_lp_draws = TRUE){
  
  #' Run one-path Pathfinder with a given list of initials in parallel
  #'
  #' @param init_ls        the list of initials 
  #' @param mc.cores       The maximum number of one-path Pathfinder to run in parallel 
  #' @param model          Stan model (compiled using rstan::stan_model())
  #' @param data           data list for call to rstan::sampling for specified model
  #' @param init_bound     the boundwith of random initials for each dimension (default = 2.0)
  #' @param N1             maxium number of iterations for L-BFGS (default = 1000)
  #' @param N_sam_DIV      number of samples used in estimating ELBO (default = 5)
  #' @param N_sam          number of samples from the approximating Normal returned
  #'                       will not do extra samples when < N_sam_DIV  (default = 5)
  #' @param factr_tol      the option factr in optim() (default = 1e2)
  #' @param lmm            the option lmm in optim() (default = 6)
  #' @param seed_list      array of random seeds for running one-path Pathfinder
  #' @param eval_lp_draws  logical; if TRUE (default) evaluate the log-densities of approximating draws
  #' @return 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  list_ind = c()
  for(i in 1:length(init_ls)){
    list_ind[[i]] <- i
  }
  out <- mclapply(list_ind, f <- function(x){
    opt_path(init = init_ls[[x]], fn = fn, gr = gr, N1 = N1,
             N_sam_DIV = N_sam_DIV, N_sam = N_sam, 
             factr_tol = factr_tol, lmm = lmm, seed = seed_list[[x]],
             eval_lp_draws = eval_lp_draws)
  }, mc.cores = mc.cores)
}

check_cond <- function(Yk, Sk){
  
  #' check whether the updates of the optimization path should be used in the 
  #' inverse Hessian estimation or not
  #' 
  #' @param Yk       update in parameters 
  #' @param Skt      update in gradient 
  #' 
  #' @return 
  
  Dk = sum(Yk * Sk)
  if(Dk == 0){
    return(FALSE)
  }else{
    thetak = sum(Yk^2) / Dk  # curvature checking
    if((Dk <= 0 | abs(thetak) > 1e12)){ #2.2*e^{-16}
      return(FALSE)
    }else{
      return(TRUE)
    }
  }
}

Form_init_Diag <- function(E0, Yk, Sk){
  
  #' Form the initial diagonal inverse Hessian in the L-BFGS update 
  #' 
  #' @param E0       initial diagonal inverse Hessian before updated
  #' @param Yk       update in parameters 
  #' @param Skt      update in gradient 
  #' 
  #' @return 
  
  Dk = sum(Yk * Sk)
  thetak = sum(Yk^2) / Dk   
  a <- (sum(E0 * Yk^2) / Dk)
  E = 1 / (a / E0 + Yk^2 / Dk - a * (Sk / E0)^2 / sum(Sk^2 / E0))
  return(E)
}

updateYS <- function(Ykt_h, Ykt, lmm){
  
  #' update the matrix for storing history of updates along optimization opath
  #' 
  #' @param Ykt_h    history of updates 
  #' @param Ykt      latest update 
  #' @param lmm      The size of history
  #' 
  #' @return
  
  if(is.null(Ykt_h)){
    Ykt_h <- matrix(Ykt, nrow = 1)
  }else if(nrow(Ykt_h) == lmm){
    Ykt_h <- rbind(Ykt_h[-1, ], Ykt)
  }else{
    Ykt_h <- rbind(Ykt_h, Ykt)
  }
  return(Ykt_h)
}

Form_N_apx_taylor <- function(x_l, g_l, Ykt_h, Skt_h, E, lmm){

  #' Returns sampling metrics of the approximating Gaussian given the history
  #' of optimization trajectory
  #'
  #' @param x_l      The point up to which the optimization path is used for approximation
  #' @param g_l      The gradient at x_l
  #' @param Ykt_h    history of updates along optimization trajectory
  #' @param Skt_h    history of updates of gradients along optimization trajectory
  #' @param E        initial diagonal inverse Hessian
  #' @param lmm      The size of history
  #'
  #' @return

  D = length(x_l)
  if(is.null(Ykt_h)){# cannot approximate Hessian
    return(NA)}

  Dk = c()
  thetak = c()
  m = nrow(Ykt_h)
  for(i in 1:m){
    Dk[i] = sum(Ykt_h[i, ] * Skt_h[i, ])
    thetak[i] = sum(Ykt_h[i, ]^2) / Dk[i]   # curvature checking
  }

  Rk = matrix(0.0, nrow = m, ncol = m)
  for(s in 1:m){
    for(i in 1:s){
      Rk[i, s] = sum(Skt_h[i, ] * Ykt_h[s, ])
    }
  }
  ninvRST = -backsolve(Rk, Skt_h)

  if( 2*m >= D){
    # directly calculate inverse Hessian and the cholesky decomposition
    Hk = diag(E, nrow = D) +
      crossprod(Ykt_h %*% diag(E, nrow = D), ninvRST) +
      crossprod(ninvRST, Ykt_h %*% diag(E, nrow = D))  +
      crossprod(ninvRST,
                (diag(Dk, nrow = m) +
                   tcrossprod(Ykt_h %*% diag(sqrt(E), nrow = D))) %*%
                  ninvRST)
    cholHk = chol(Hk)
    logdetcholHk = determinant(cholHk)$modulus

    x_center = c(x_l - Hk %*% g_l) # consider the first term in normal approximation

    sample_pkg <- list(label = "full", cholHk = cholHk,
                       logdetcholHk = logdetcholHk,
                       x_center = x_center)

  } else {
    # use equation ?? to sample
    Wkbart = rbind(Ykt_h %*% diag(sqrt(E)),
                   ninvRST %*% diag(sqrt(1 / E)))
    Mkbar = rbind(cbind(matrix(0.0, nrow = m, ncol = m), diag(m)),
                  cbind(diag(m),
                        (diag(Dk, nrow = m) +
                           tcrossprod(Ykt_h %*% diag(sqrt(E))))))
    qrW = qr(t(Wkbart))
    Qk = qr.Q(qrW)
    Rkbar = qr.R(qrW)
    Rktilde = chol(Rkbar %*% Mkbar %*% t(Rkbar) + diag(nrow(Rkbar)))
    logdetcholHk = sum(log(diag(Rktilde))) + 0.5 * sum(log(E))

    ninvRSTg = ninvRST %*% g_l
    x_center = c(x_l -   # consider the first term in taylor expansion
      (E*g_l + E*crossprod(Ykt_h, ninvRSTg) +
      crossprod(ninvRST, Ykt_h %*% (E*g_l))  +
      crossprod(ninvRST,
                (diag(Dk, nrow = m) +
                   tcrossprod(Ykt_h %*% diag(sqrt(E), nrow = D))) %*%
                  ninvRSTg)))

    sample_pkg <- list(label = "sparse", theta_D = 1 / E,
                       Qk = Qk, Rktilde = Rktilde,
                       logdetcholHk = logdetcholHk,
                       Mkbar = Mkbar,
                       Wkbart = Wkbart,
                       x_center = x_center)
  }
  return(sample_pkg)
}

est_DIV <- function(sample_pkg, N_sam, fn, label = "ELBO"){
  
  #' estimate divergence based on Monte Carlo samples given the output of 
  #' subroutine Form_N_apx
  #' 
  #' @param sample_pkg  The output of subroutine Form_N_apx
  #' @param N_sam       number of samples used in estimating divergence
  #' @param fn          negative of log density
  #' @param label       type of divergence, (default = "ELBO")
  #' 
  #' @return 
  
  D <- length(sample_pkg$x_center)
  fn_draws <-  rep(Inf, N_sam)
  lp_approx_draws <- rep(0.0, N_sam)
  repeat_draws <- matrix(0, nrow = D, ncol = N_sam)
  draw_ind = 1
  fn_calls_DIV = 0
  f_test_DIV = -Inf
  
  for(l in 1:(2*N_sam)){
    if(sample_pkg$label == "full"){
      u = rnorm(D)
      u2 = crossprod(sample_pkg$cholHk, u) + c(sample_pkg$x_center)
    }else{
      u = rnorm(D)
      u1 = crossprod(sample_pkg$Qk, u)
      u2 = diag(sqrt(1 / sample_pkg$theta_D)) %*%
        (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, u1) + 
           (u - sample_pkg$Qk %*% u1)) + c(sample_pkg$x_center)
    }
    # skip bad samples
    skip_flag = FALSE
    tryCatch(f_test_DIV <- fn(u2),
             error = function(e) { skip_flag <<- TRUE})
    if(is.nan(f_test_DIV)){skip_flag <<- TRUE}
    if(skip_flag){
      next
    } else {
      fn_draws[draw_ind] <- f_test_DIV
      lp_approx_draws[draw_ind] <- - sample_pkg$logdetcholHk - 0.5 * sum(u^2) - 
        0.5 * D * log(2 * pi)
      repeat_draws[, draw_ind] <- u2
      draw_ind = draw_ind + 1
    }
    fn_calls_DIV = fn_calls_DIV + 1
    if(draw_ind == N_sam + 1){break}
  }
  
  ### Divergence estimation ###
  ELBO <- -mean(fn_draws) - mean(lp_approx_draws)
  if(is.nan(ELBO)){ELBO <- -Inf}
  if(label == "ELBO"){
    DIV <- ELBO
  }else if(label == "lIKL"){
    DIV <- ELBO + log(mean(exp(-fn_draws - lp_approx_draws - ELBO)))    # log Inclusive-KL E[p(x_i)/q(x_k)]
  }else if(label == "lADIV"){
    DIV <- 0.5 * ELBO + 
      log(mean(exp(0.5 * (-fn_draws - lp_approx_draws - ELBO)))) # log alpha-divergence E[(p(x_i)/q(x_k))^alpha], e.g. with 1/2
  }else if(label == "lCDIV"){
    DIV <- 2.0 * ELBO + 
      log(mean(exp(2.0 * (-fn_draws - lp_approx_draws - ELBO))))# log Chi^2-divergence is alpha-divergence with alpha=2
  } else{
    stop("The divergence is misspecified")
  }
  
  if(is.nan(DIV)){DIV <- -Inf}
  if(is.infinite(DIV)){DIV <- -Inf}
  
  if(any(is.nan(fn_draws))){
    fn_draws[is.nan(fn_draws)] <- Inf
  }
  
  return(list(DIV = DIV, 
              repeat_draws = repeat_draws, fn_draws = fn_draws, 
              lp_approx_draws = lp_approx_draws, fn_calls_DIV = fn_calls_DIV))
}

Sam_N_apx <- function(sample_pkg, N_sam){
  
  #' Generate N_sam samples from the approximating Normal given the output of 
  #' subroutine Form_N_apx
  #' 
  #' @param sample_pkg The output of subroutine Form_N_apx
  #' @param N_sam      Number of samples 
  #' 
  #' @return \item{samples}{The samples from the approximating Normal, 
  #'   each column stores a sample} 
  #'   \item{lp_approx_draws}{The log-density of generated samples under 
  #'   approximating Normal} 
  
  
  D <- length(sample_pkg$x_center)
  lp_approx_draws <- rep(0.0, N_sam)
  repeat_draws <- matrix(0, nrow = D, ncol = N_sam)
  
  if(sample_pkg$label == "full"){
    u = matrix(rnorm(D * N_sam), nrow = D)
    u2 = crossprod(sample_pkg$cholHk, u) + c(sample_pkg$x_center)
  }else{
    u = matrix(rnorm(D * N_sam), nrow = D)
    u1 = crossprod(sample_pkg$Qk, u)
    u2 = diag(sqrt(1 / sample_pkg$theta_D)) %*%
      (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, u1) + 
         (u - sample_pkg$Qk %*% u1)) + c(sample_pkg$x_center)
  }
  
  lp_apx_draws <- - sample_pkg$logdetcholHk - 0.5 * colSums(u^2) - 
    0.5 * D * log(2*pi)
  
  return(list(samples = u2, lp_apx_draws = lp_apx_draws))
}


params_only <- function(path) {
  
  #' Return optimization path with last column (objective function value)
  #' removed.
  #'
  #' @param path   optimization path with last column for objective
  #'               function value
  #' @return
  
  N <- dim(path)[1]
  D <- dim(path)[2]
  path[1:N, 1:(D - 1)]
}

fit_info_DIV <- function(param_path) {
  # extract fitting information #
  DIV = param_path$DIV_save$DIV
  
  return(DIV)
}


Imp_Resam_WOR <- function(param_path, n_inits, seed = 123){
  
  #' PSIS-IR without replacement for multi-path Pathfinder 
  #' Return n_inits distinct samples
  #'  
  #' @param param_path output of function opt_path_stan_parallel()
  #' @param n_inits    number of distinct samples
  #' @param seed       random seed for importance resampling
  
  # remove one-path Pathfinder that returns mode
  check <- sapply(param_path, f <- function(x){
    work <- (x$status == "samples")
    work
  })
  
  filter_mode <- c(1:length(check))[check]
  
  param_path <- param_path[filter_mode]
  
  ## extract samples and log ratios ##
  samples <- do.call("cbind", lapply(param_path, extract_samples))
  lrms <- c(sapply(param_path, extract_log_ratio))
  
  ## take off samples with infinite log ratios
  finit_ind <- is.finite(lrms)
  samples <- samples[, finit_ind]
  lrms <- lrms[finit_ind]
  
  ## compute the importance weight ##
  sample_weights_psis <- suppressWarnings(weights(psis(lrms, r_eff = NA),
                                                  log = FALSE))
  #sample_weights_IS <- exp(lrms - max(lrms))/sum(exp(lrms - max(lrms)))
  
  ## Importance resampling ##
  set.seed(seed)
  sample_ind <-sample.int(ncol(samples), size = n_inits, replace = FALSE, 
                          prob = sample_weights_psis)
  
  return(samples[, sample_ind])
}

Imp_Resam_WR <- function(param_path, n_sam, seed = 123){
  
  #' PSIS-IR with replacement for multi-path Pathfinder 
  #' Return n_sam samples
  #'  
  #' @param param_path output of function opt_path_stan_parallel()
  #' @param n_sam    number of distinct samples
  #' @param seed       random seed for importance resampling
  
  # remove one-path Pathfinder that returns mode
  check <- sapply(param_path, f <- function(x){
    work <- (x$status == "samples")
    work
  })
  
  filter_mode <- c(1:length(check))[check]
  
  param_path <- param_path[filter_mode]
  
  ## extract samples and log ratios ##
  samples <- do.call("cbind", lapply(param_path, extract_samples))
  lrms <- c(sapply(param_path, extract_log_ratio))
  
  ## take off samples with infinite log ratios
  finit_ind <- is.finite(lrms)
  samples <- samples[, finit_ind]
  lrms <- lrms[finit_ind]
  
  ## compute the importance weight ##
  sample_weights_psis <- suppressWarnings(weights(psis(lrms, r_eff = NA),
                                                  log = FALSE))
  #sample_weights_IS <- exp(lrms - max(lrms))/sum(exp(lrms - max(lrms)))
  
  ## Importance resampling ##
  set.seed(seed)
  sample_ind <-sample.int(ncol(samples), size = n_sam, replace = TRUE, 
                          prob = sample_weights_psis)
  
  return(samples[, sample_ind])
}

Imp_Resam_Each <- function(param_path, seed){
  
  #' PSIS for one-path Pathfinder 
  #'  
  #' @param param_path output of function opt_path_stan_parallel()
  #' @param seed       random seed for importance resampling
  
  set.seed(seed)
  
  # sample using importance sampling
  EX_ests <- lapply(param_path, 
                    f <- function(x){
                      if(x$status == "mode"){ # return mode
                        EX_est <- x$y[nrow(x$y), -ncol(x$y)]
                      } else {
                        lrms <- extract_log_ratio(x)
                        samples <- extract_samples(x)
                        finit_ind <- is.finite(lrms) 
                        samples <- samples[, finit_ind]
                        lrms <- lrms[finit_ind]
                        
                        psis_w <- suppressWarnings(
                          weights(psis(lrms, r_eff = NA), log = FALSE))
                        pick_ind <- sample(1:length(psis_w),
                                           replace = TRUE,
                                           size = 1, prob = psis_w)
                        EX_est <- samples[, pick_ind]
                      }
                      return(EX_est)})
  Sample_Each = do.call("cbind", EX_ests)
  return(Sample_Each)
}

random_sample_Each <- function(param_path, seed){
  
  #' random sample for one-path Pathfinder 
  #'  
  #' @param param_path output of function opt_path_stan_parallel()
  #' @param seed       random seed 
  
  set.seed(seed)
  
  # sample using importance sampling
  EX_ests <- lapply(param_path, 
                    f <- function(x){
                      if(x$status == "mode"){ # return mode
                        EX_est <- x$y[nrow(x$y), -ncol(x$y)]
                      } else {
                        #lrms <- extract_log_ratio(x)
                        samples <- extract_samples(x)
                        #finit_ind <- which(is.finite(lrms))
                        #pick_ind <- sample(finit_ind, size = 1)
                        pick_ind <- sample.int(ncol(samples), size = 1)
                        EX_est <- samples[, pick_ind]
                      }
                      return(EX_est)})
  Sample_Each = do.call("cbind", EX_ests)
  return(Sample_Each)
}


extract_samples <- function(param_path){
  sample <- param_path$DIV_save$repeat_draws
  return(sample)
}

extract_samples_center <- function(param_path){
  sample <- param_path$sample_pkg_save$x_center
  return(sample)
}

extract_lps <- function(param_path){
  lps <- -param_path$DIV_save$fn_draws
  return(lps)
}

constrain_draws <- function(unconstrain_draws, posterior){
  lapply(1:ncol(unconstrain_draws), 
         f <- function(s){
           constrain_pars(posterior, unconstrain_draws[, s])
         })
}

extract_log_ratio <- function(param_path){
  if(length(param_path$DIV_save$fn_draws)!=
     length(param_path$DIV_save$lp_approx_draws)){
    stop('length of log-densities for approximate draws does not match with the sample size, set option eval_lp_draws at TRUE before running importance resampling')
  }
  lrm <- - param_path$DIV_save$fn_draws -
    param_path$DIV_save$lp_approx_draws
  return(lrm)
}


filter_samples_resam <- function(param_path, n_inits, seed = 123){
  
  #' SIR with mixture Gaussian 
  
  # remove the failed Pathfinder
  check <- sapply(param_path, f <- function(x){
    work <- TRUE
    tryCatch(
      lps <- extract_lps(x), error = function(e) { work <<- FALSE})
    work
  })
  filter_mode <- c(1:length(check))[check]
  
  param_path <- param_path[filter_mode]
  
  lps <- c(sapply(param_path, extract_lps))
  samples <- lapply(param_path, extract_samples)
  
  J <- length(param_path)
  ns <- ncol(samples[[1]])
  lp_approx_M <- matrix(0, nrow = J*ns, ncol = J) # save log-densities 
  
  for(indI in 1:J){  #sample
    for(indJ in 1:J){ #pathfinder
      if (indJ == indI){
        lp_approx_M[((indI-1)*ns + 1):(indI * ns), indJ] <- 
          param_path[[indJ]]$DIV_save$lp_approx_draws
      } else {
        lps_tem <- rep(0, ns)
        if(param_path[[indJ]]$sample_pkg_save$label == "full"){
          u <- forwardsolve(param_path[[indJ]]$sample_pkg_save$cholHk,
                            (samples[[indI]] -  
                               param_path[[indJ]]$sample_pkg_save$x_center),
                            transpose = TRUE, upper.tri = TRUE)
          lps_tem = - param_path[[indJ]]$sample_pkg_save$logdetcholHk -
            0.5 * colSums(u^2)
        }else{
          u1 <- (samples[[indI]] - param_path[[indJ]]$sample_pkg_save$x_center) * 
            sqrt(param_path[[indJ]]$sample_pkg_save$theta_D)
          u2 = crossprod(param_path[[indJ]]$sample_pkg_save$Qk, u1)
          norm_u <- 
            colSums(forwardsolve(param_path[[indJ]]$sample_pkg_save$Rktilde, u2,
                                 transpose = TRUE, upper.tri = TRUE)^2) + 
            colSums(u1^2) - colSums(u2^2)
          lps_tem = - param_path[[indJ]]$sample_pkg_save$logdetcholHk -
            0.5 * norm_u
        }
        lp_approx_M[((indI-1)*ns + 1):(indI * ns), indJ] <- lps_tem
      }
    }
  }
  
  lp_approx <- apply(lp_approx_M, 1, logSumExp)
  lrms <- lps - lp_approx
  
  sample_weights <- suppressWarnings(weights(psis(lrms, r_eff = NA), 
                                             log = FALSE))
  
  samples_M <- matrix(unlist(samples), nrow = nrow(samples[[1]]))
  
  set.seed(seed)
  sample_ind <- sample(1:length(sample_weights), replace = TRUE,
                       size = n_inits, prob = sample_weights)
  print(table(sample_ind))
  return(samples_M[, sample_ind])
}


get_opt_tr <- function(opath){
  
  ###
  #' function for retreveing optimization trajectories
  #' 
  
  lp_ind = ncol(opath[[1]]$y)
  opt_tr <- c()
  ind_tr <- c()
  tr_id <- c()
  for(l in 1:length(opath)){
    opt_tr = rbind(opt_tr, opath[[l]]$y[, 1:(lp_ind - 1)])
    ind_tr = c(ind_tr, 
               1:nrow(opath[[l]]$y[, 1:(lp_ind - 1)]))
    tr_id = c(tr_id, 
              rep(l, nrow(opath[[l]]$y[, 1:(lp_ind - 1)])))
    
  }
  return(list(opt_tr = opt_tr, ind_tr = ind_tr, tr_id = tr_id))
}


Form_N_apx <- function(x_center, Ykt_h, Skt_h, E, lmm){

  #' Returns sampling metrics of the approximating Gaussian given the history
  #' of optimization trajectory
  #'
  #' @param x_center The center of the approximating Gaussian
  #' @param Ykt_h    history of updates along optimization trajectory
  #' @param Skt_h    history of updates of gradients along optimization trajectory
  #' @param E        initial diagonal inverse Hessian
  #' @param lmm      The size of history
  #'
  #' @return

  D = length(x_center)
  if(is.null(Ykt_h)){# cannot approximate Hessian
    return(NA)}

  Dk = c()
  thetak = c()
  m = nrow(Ykt_h)
  for(i in 1:m){
    Dk[i] = sum(Ykt_h[i, ] * Skt_h[i, ])
    thetak[i] = sum(Ykt_h[i, ]^2) / Dk[i]   # curvature checking
  }

  Rk = matrix(0.0, nrow = m, ncol = m)
  for(s in 1:m){
    for(i in 1:s){
      Rk[i, s] = sum(Skt_h[i, ] * Ykt_h[s, ])
    }
  }
  ninvRST = -backsolve(Rk, Skt_h)

  if( 2*m >= D){
    # directly calculate inverse Hessian and the cholesky decomposition
    Hk = diag(E, nrow = D) +
      crossprod(Ykt_h %*% diag(E, nrow = D), ninvRST) +
      crossprod(ninvRST, Ykt_h %*% diag(E, nrow = D))  +
      crossprod(ninvRST,
                (diag(Dk, nrow = m) +
                   tcrossprod(Ykt_h %*% diag(sqrt(E), nrow = D))) %*%
                  ninvRST)
    cholHk = chol(Hk)
    logdetcholHk = determinant(cholHk)$modulus

    sample_pkg <- list(label = "full", cholHk = cholHk,
                       logdetcholHk = logdetcholHk,
                       x_center = x_center)

  } else {
    # use equation ?? to sample
    Wkbart = rbind(Ykt_h %*% diag(sqrt(E)),
                   ninvRST %*% diag(sqrt(1 / E)))
    Mkbar = rbind(cbind(matrix(0.0, nrow = m, ncol = m), diag(m)),
                  cbind(diag(m),
                        (diag(Dk, nrow = m) +
                           tcrossprod(Ykt_h %*% diag(sqrt(E))))))
    qrW = qr(t(Wkbart))
    Qk = qr.Q(qrW)
    Rkbar = qr.R(qrW)
    Rktilde = chol(Rkbar %*% Mkbar %*% t(Rkbar) + diag(nrow(Rkbar)))
    logdetcholHk = sum(log(diag(Rktilde))) + 0.5 * sum(log(E))

    sample_pkg <- list(label = "sparse", theta_D = 1 / E,
                       Qk = Qk, Rktilde = Rktilde,
                       logdetcholHk = logdetcholHk,
                       Mkbar = Mkbar,
                       Wkbart = Wkbart,
                       x_center = x_center) #inv_metric_est is returned for estimating mass matrix
  }
  return(sample_pkg)
}

