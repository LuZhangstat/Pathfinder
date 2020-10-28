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
# For N iterations, the algorithm makes N calls to the base
# R function optim(), with max iteration set to n for each n in 1:N.
# This allows programmatic access to the iterations of optim(), but is
# O(N^2). Thanks to Calvin Whealton for the approach:
#   https://stackoverflow.com/a/46904780).
#
# @param x_init initial parameter values
# @param fn objective function
# @param gr gradient function of objective
# @param N total number of iterations (default 25)
# @return optimization path, matrix with N rows, each of which is a
#         point on the optimization path followed by objective


opt_path <- function(init, fn, gr, lp_f, 
                     init_bound = 2.0,
                     N1 = 100, N_mode_max = 20, 
                     N_sam = 100) {
  D <- length(init)
  #init <- runif(D, -init_bound, init_bound)
  c = 1.0; rho = 1.0
  # preallocation
  E_lp <- c()
  E <- c()
  #ill_cond = c()    #0.0 no problem, 1.0 not pd Hessian
  step_count = c()
  
  y <- matrix(NA, nrow = 1, ncol = D + 2)
  y[1, 1:D] <- init
  y[1, D + 1] <- fn(init)
  y[1, D + 2] <- lp_f(init)
  n0 = 1
  abs_tol = 1e-2; accept_tol = 1e-13;
  factor_tol = 1e-5;
  
  for (j in 1:N_mode_max){
    cat(j, "\t")
    if(j == 1){
      RHR_fail <- FALSE
      tryCatch(RHR_fit <- RHR(fn, gr, y[n0, 1:D], sigma0 = 1.0, max_iter = N1, 
                              abs_tol = abs_tol, factor_tol = factor_tol,
                              accept_tol = accept_tol), 
               error = function(e) { RHR_fail <<- TRUE})
      if(RHR_fail){
        for(l in 1:40){
          print("\n reinitialize \n")
          y[n0, 1:D] <- runif(D, -init_bound, init_bound)
          y[n0, D + 1] <- fn(y[n0, 1:D])
          y[n0, D + 2] <- lp_f(y[n0, 1:D])
          RHR_fail <- FALSE
          tryCatch(RHR_fit <- RHR(fn, gr, y[n0, 1:D], sigma0 = 1.0, max_iter = N1, 
                                  abs_tol = abs_tol, factor_tol = factor_tol,
                                  accept_tol = accept_tol), 
                   error = function(e) { RHR_fail <<- TRUE})
          if(!RHR_fail){break}
        }
      }
    } else { #j > 1
      RHR_fail <- FALSE
      tryCatch(RHR_fit <- RHR(fn, gr, y[n0, 1:D], sigma0 = 1.0, max_iter = N1, 
                              abs_tol = abs_tol, factor_tol = factor_tol,
                              accept_tol = accept_tol), 
               error = function(e) { RHR_fail <<- TRUE})
    }
    
    
    if(RHR_fail){
      if(j == 1){
        E_lp[1] <- NA
        E[1] <- NA
        y = y
        return(list(E_lp = E_lp, E = E, y = y))
      }
      break} # fail in j > 1, break
    
    lp_ls <- sapply(RHR_fit$x_list, lp_f)
    x_track <- matrix(unlist(RHR_fit$x_list), byrow = TRUE, 
                      nrow = length(RHR_fit$x_list), ncol = D)
    y = rbind(y, 
              cbind(x_track, RHR_fit$fn_ls, lp_ls))
    
    # save results
    step_count[j] <- sum(RHR_fit$steps)
    
    
    # simulation samples of approximated posterior distribution and estimate E_lp
    
    mode_samples <- 
      RHR_monte_carlo(lp_f, 
                      fn, 
                      D, 
                      rk = RHR_fit$rank, 
                      Rk = RHR_fit$Rk, 
                      Zk = RHR_fit$Zk, 
                      xk = RHR_fit$x_list[[length(RHR_fit$x_list)]],
                      sigmak = RHR_fit$sigmak, 
                      N_sam)
    
    E_lp[j] = mode_samples$E_lp
    E[j] = mode_samples$E
    
    if(mode_samples$next_flag && (j < N_mode_max)){
      n0 <- nrow(y) + 1
      y = rbind(y, c(mode_samples$x0_up, fn(mode_samples$x0_up), 
                     lp_f(mode_samples$x0_up)))
    } else {
      break
    }
  }
  #E_lp
  return(list(E_lp = E_lp, E = E, y = y, step_count = step_count))
}


# experimental code
# tt <- capture.output(z <- optim(par = x_init,
#                           fn = function(x) -fn(x),  # negate for maximization
#                           gr = function(x) -gr(x),
#                           method = "L-BFGS-B",
#                           control = list(maxit = N, factr = 1e8,
#                                          trace = 6, REPORT = 1)))
# 
# tt2 <- sapply(strsplit(tt[seq(23, (length(tt) -14), by = 4)], split = " "), 
#        f <- function(x){as.numeric(x[c(-1, -2)])})
  

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
opt_path_stan <- function(model, data, N1, N_mode_max, N_sam, init_bound = 2) {
  # require chains > 0, iter > 0, so use Fixed_param to avoid work
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, adjust_transform = FALSE)
  lp_f <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE, 
                                 gradient = TRUE)[1]
  out <- opt_path(init, fn, gr, lp_f, N1, N_mode_max, N_sam)
  return(out)
}
#opt_path_stan(model, data, N1, N_rep, init_bound = 2)

opt_path_stan_parallel <- function(seed_list, mc.cores, model, data, 
                                   N1, N_mode_max, N_sam, init_bound){
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, adjust_transform = FALSE)
  lp_f <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE, 
                                   gradient = TRUE)[1]
  MC = length(seed_list)
  init = c()
  for(i in 1:MC){
    set.seed(seed_list[i])
    init[[i]] <- runif(D, -init_bound, init_bound)
  }
  out <- mclapply(init, opt_path, fn = fn, gr = gr, lp_f = lp_f, N1 = N1, 
                  N_mode_max = N_mode_max, N_sam = N_sam, mc.cores = mc.cores)
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
  mode_ind = length(param_path$E)
  E_target = param_path$E[mode_ind]
  mark_index = which.min(abs(param_path$y[, (ncol(param_path$y) - 1)] - E_target))
  mark_index2 = which.min(abs(param_path$y[, (ncol(param_path$y))] - 
                                param_path$E_lp[mode_ind]))
  
  return(mark_index2)
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

