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
opt_path <- function(x_init, fn, gr, N = 25, lp_f) {
  D <- length(x_init)
  y <- matrix(NA, nrow = N + 1, ncol = D + 1)
  y[1, 1:D] <- x_init
  y[1, D + 1] <- fn(x_init)
  for (n in 1:N) {
    break_opt <- FALSE
    tryCatch(z <- optim(par = x_init,
                        fn = function(x) -fn(x),  # negate for maximization
                        gr = function(x) -gr(x),
                        method = "L-BFGS-B",
                        control = list(maxit = n#, factr = 1e10#, ndeps = 1e-8 #, 
                                       #trace = 6, REPORT = 1 
                                       )), 
             error = function(e) { break_opt <<- TRUE})
    if(break_opt) { 
      print("Error in obtaining optimization path.")
      return(y[1:n, 1:(D + 1)])
    }
    
    y[n + 1, 1:D] <- z$par
    y[n + 1, D + 1] <- lp_f(z$par)
    # break if no change in objective
    printf("n = %4d, last = %f;   this = %f",
           n, y[n, D + 1], y[n + 1, D + 1])
    print(z$counts)
    if (y[n, D + 1] == y[n + 1, D + 1]) {
      return(y[1:n, 1:(D + 1)])
    }
  }
  y
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
opt_path_stan <- function(model, data, N = 25, init_bound = 2) {
  # require chains > 0, iter > 0, so use Fixed_param to avoid work
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) grad_log_prob(posterior, theta, adjust_transform = FALSE)
  lp_f <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE, 
                                 gradient = TRUE)[1]
  out <- opt_path(init, fn, gr, N, lp_f)
  return(out)
}

opt_path_stan_parallel <- function(seed_list, mc.cores, 
                                   model, data, N, init_bound){
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
  out <- mclapply(init, opt_path, fn = fn, gr = gr, N = N, lp_f = lp_f,
                  mc.cores = mc.cores)
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

lp_draws <- function(posterior, init_param_unc, int_time, fn, gr, stepsize) {
  
  theta_up <- init_param_unc
  counts_inc = 0
  trapped = 0
  mark = 0
  draws <- c()
  g0 <- gr(theta_up)
  lp0 <- fn(theta_up)
  rho0 = (sqrt(length(theta_up)) / sqrt(sum(g0^2)))*g0
  rho <- rho0 + 0.5 * stepsize * g0
  break_flag <- FALSE
  for (l in 1:int_time){
    theta_up <- theta_up + stepsize * rho
    tryCatch(draws[l] <- fn(theta_up),
             error = function(e) { break_flag <<- TRUE })
    if(break_flag | is.na(draws[l])){break}
    if(draws[l] > lp0){ counts_inc = counts_inc + 1 } # record increase lp or not
    if(draws[l] < lp0 - 0.4 * sum(rho0^2) && 
       draws[l] > lp0 - 0.6 * sum(rho0^2)){ trapped = trapped + 1 } # record trapped or not
    if(l < int_time){
      tryCatch(g1 <- gr(theta_up),
              error = function(e) { break })
      if(break_flag | any(is.na(g1))){break}
      rho <- rho + stepsize * g1
    }
  }
  if(counts_inc > 2 && trapped > 0 && 
     (max(draws, na.rm = TRUE) -lp0) < sum(rho0^2)){ 
    mark = 1 # record mark or not, mark the points that reach from the bottom of a curvature
  } 
  return(list(mark = mark, draws = c(lp0, draws), counts_inc = counts_inc,
              trapped = trapped))
}

lp_draws_2 <- function(posterior, init_param_unc, int_time, fn, gr, stepsize) {
  
  theta_up <- init_param_unc
  counts_inc = 0
  trapped = 0
  mark = 0
  draws <- c()
  g0 <- gr(theta_up)
  lp0 <- fn(theta_up)
  rho0 = -(sqrt(length(theta_up)) / sqrt(sum(g0^2)))*g0
  rho <- rho0 + 0.5 * stepsize * g0
  break_flag <- FALSE
  for (l in 1:int_time){
    theta_up <- theta_up + stepsize * rho
    tryCatch(draws[l] <- fn(theta_up),
             error = function(e) { break_flag <<- TRUE })
    if(break_flag | is.na(draws[l])){break}
    if(draws[l] > lp0){ counts_inc = counts_inc + 1 } # record increase lp or not
    if(draws[l] < lp0 - 0.4 * sum(rho0^2) && 
       draws[l] > lp0 - 0.6 * sum(rho0^2)){ trapped = trapped + 1 } # record trapped or not
    if(l < int_time){
      tryCatch(g1 <- gr(theta_up),
               error = function(e) { break })
      if(break_flag | any(is.na(g1))){break}
      rho <- rho + stepsize * g1
    }
  }
  if(counts_inc > 2 && trapped > 0 && 
     (max(draws, na.rm = TRUE) -lp0) < 0.5 * sum(rho0^2)){ 
    mark = 1 # record mark or not, if mark then break the dynamic trajectories
    #break
  } 
  return(list(mark = mark, draws = draws, counts_inc = counts_inc,
              trapped = trapped))
}

is_typical <- function(model, data, init_param_unc, M, int_time, lp_0) {
  posterior <- to_posterior(model, data)
  init_fun <- function(chain_id) constrain_pars(posterior, init_param_unc)
  # find a stepsize #
  fit_0 <- sampling(model, data = data, init = init_fun,
                    chains = 1, iter = 2, warmup = 1, refresh = 0,
                    control = list(metric = "unit_e",
                                   adapt_engaged = TRUE,
                                   max_treedepth = 1 #,
                                   # <-  stepsize
                    ),
                    save_warmup = TRUE)
  stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"]/2
  
  # generate hamiltonian dynamics
  fn <- function(theta) log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) grad_log_prob(posterior, theta, 
                                      adjust_transform = FALSE)
  lps <- lp_draws(posterior, init_param_unc, int_time, fn, gr, stepsize)
  #lps2 <- lp_draws_2(posterior, init_param_unc, int_time, fn, gr, stepsize)
  
  return(list(lps = lps))
}

find_typical <- function(param_path, model, data, M = 4, int_time = 30) {
  typical_index <- c()          # return the index of sample in param_path that is identified as a good initial
  N <- dim(param_path)[1]
  D <- dim(param_path)[2] - 1   # includes objective in last position
  record_lp_draws <- list()
  for (n in 1:N) {
    record_lp_draws[[n]] <- list()
    typl_res <- is_typical(model, data, param_path[n, 1:D], M, int_time,
                                  param_path[n, D + 1])
    record_lp_draws[[n]] <- typl_res$lps$draws
    printf("n = %3d  counts_inc1 = %3d, trapped1 = %3d",
           n, typl_res$lps$counts_inc, typl_res$lps$trapped)
    
    if (typl_res$lps$mark == 1){
      print(param_path[n, 1:(D + 1)], digits = 2)
      typical_index <- c(typical_index, n)
    }
  }
  return(list(typical_index = typical_index, 
              record_lp_draws = record_lp_draws))
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

