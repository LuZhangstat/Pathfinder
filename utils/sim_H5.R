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


opt_path <- function(init, fn, gr, lp_f, N1 = 30, N_mode_max = 20, N_sam = 100) {
  D <- length(init)
  #init <- runif(D, -init_bound, init_bound)
  c = 1.0; rho = 1.0
  # preallocation
  LBFGS_H_list <- list()
  E <- c()
  E_lp <- c()
  Mass <- c()
  ill_cond = c()    #0.0 no problem, 1.0 not pd Hessian
  fn_count = c()
  gr_count = c()
  
  D <- length(init)
  y <- matrix(NA, nrow = N1*N_mode_max + 1, ncol = D + 2)
  y[1, 1:D] <- init
  y[1, D + 1] <- fn(init)
  y[1, D + 2] <- lp_f(init)
  n0 = 1
  for (j in 1:N_mode_max){
    for(n in 1:N1){
      break_opt <- FALSE
      tryCatch(z <- optim(par = y[n0, 1:D],
                          fn = function(x) -fn(x),  # negate for maximization
                          gr = function(x) -gr(x),
                          method = "L-BFGS-B",
                          hessian = TRUE,
                          control = list(maxit = n, factr = 1e9#, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        if(n == 1){
          print("Error in obtaining optimization path.")
          return(list(LBFGS_H_list = LBFGS_H_list,
                      E = E, E_lp = E_lp, Mass = Mass,
                      ill_cond = ill_cond, y = y[1:(n0), ]))
        }
        print("Error in obtaining optimization path.")
        #return(y[1:(n0 + n), 1:(D + 1)])
        n = n - 1
        break
      }
      y[n + n0, 1:D] <- z$par
      y[n + n0, D + 1] <- -z$value
      y[n + n0, D + 2] <- lp_f(z$par)
      
      # break if no change in objective
      printf("n = %4d, last = %f;   this = %f",
             n , y[n + n0 - 1, D + 2], y[n + n0, D + 2])
      print(z$counts)
      if (y[n + n0 - 1, D + 1] == y[n + n0, D + 1]) {
        n = n - 1; break
      }
    }
    # save results
    LBFGS_H_list[[j]] <- z$hessian
    fn_count[j] <- z$counts[1]
    gr_count[j] <- z$counts[2]
    
    # simulation samples of approximated posterior distribution
    ill_condition <- FALSE
    ill_cond[j] = 0.0
    tryCatch(cholH <- chol(z$hessian), 
             error = function(e) { ill_condition <<- TRUE})
    if(ill_condition){
      ill_cond[j] = 1.0
      return(list(LBFGS_H_list = LBFGS_H_list,
                  E = E, E_lp = E_lp, Mass = Mass,
                  ill_cond = ill_cond, y = y[1:(n0 + n), ]))
    }
    u = matrix(rnorm(D * N_sam), ncol = N_sam)
    u2 = solve(cholH)%*% u 
    fn_draws <-  c()
    ind_draws <- c()
    lp_f_draws <- c()
    for(l in 1:ncol(u)){
      skip_flag = FALSE
      tryCatch(f_test <- fn(y[n + n0, 1:D] + c * rho^j * u2[, l]),
               error = function(e) { skip_flag <<- TRUE})
      if(skip_flag){next}
      else {
        ind_draws <- c(ind_draws, l)
        fn_draws <- c(fn_draws, f_test)
        lp_f_draws <- c(lp_f_draws, lp_f(y[n + n0, 1:D] + c * rho^j * u2[, l]))
      }
    }
    E[j] = mean(fn_draws)
    E_lp[j] = mean(lp_f_draws)
    Mass[j] = 0.0
    next_flag = FALSE
    if(max(fn_draws) > y[n+n0, D+1]){next_flag <- TRUE}
    #cat(next_flag)
    if(next_flag && (j < N_mode_max)){
      init <- y[n0 + n, 1:D] + c * rho^j * u2[, ind_draws[which.max(fn_draws)]]
      n0 <- n0 + n
    } else {
      return(list(LBFGS_H_list = LBFGS_H_list,
                  E = E, E_lp = E_lp, Mass = Mass,
                  ill_cond = ill_cond, y = y[1:(n0 + n), ]))
    }
  }
  return(list(LBFGS_H_list = LBFGS_H_list,
              E = E, E_lp = E_lp, Mass = Mass,
              ill_cond = ill_cond, y = y[1:(n0 + n), ]))
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
  fn <- function(theta) log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) grad_log_prob(posterior, theta, adjust_transform = FALSE)
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
  fn <- function(theta) log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) grad_log_prob(posterior, theta, adjust_transform = FALSE)
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
  mark_index = which.min(abs(param_path$y[, ncol(param_path$y) - 1] - E_target))
  return(mark_index)
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

