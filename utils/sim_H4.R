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
# opt_path <- function(x_init, fn, gr, N1 = 30, N_rep = 6, lp_f) {
#   D <- length(x_init)
#   y <- matrix(NA, nrow = N1 + N2 + 1, ncol = D + 1)
#   y[1, 1:D] <- x_init
#   y[1, D + 1] <- lp_f(x_init)
#   for (n in 1:N1) {
#     break_opt <- FALSE
#     tryCatch(z <- optim(par = x_init,
#                         fn = function(x) -fn(x),  # negate for maximization
#                         gr = function(x) -gr(x),
#                         method = "L-BFGS-B",
#                         hessian = TRUE,
#                         control = list(maxit = n, factr = 1e10#, ndeps = 1e-8 #, 
#                                        #trace = 6, REPORT = 1 
#                                        )), 
#              error = function(e) { break_opt <<- TRUE})
#     if(break_opt) { 
#       print("Error in obtaining optimization path.")
#       return(y[1:n, 1:(D + 1)])
#     }
#     
#     y[n + 1, 1:D] <- z$par
#     y[n + 1, D + 1] <- lp_f(z$par)
#     
#     # break if no change in objective
#     printf("n = %4d, last = %f;   this = %f",
#            n, y[n, D + 1], y[n + 1, D + 1])
#     print(z$counts)
#     if (y[n, D + 1] == y[n + 1, D + 1]) {
#       return(y[1:n, 1:(D + 1)])
#     }
#   }
#   sk = (1 / sqrt(sum(gr(x_init)^2)))/100
#   #sk = 
#   vk = rep(0.0, D)
#   vk_1 = vk
#   mk = 1.0
#   fn_old = fn(z$par)
#   for (n in (N1 + 1):(N1 + N2)){
#     vk_1 = y[n, 1:D] + sk * gr(y[n, 1:D])
#     y[n + 1, 1:D] = vk_1 + (mk - 1.0) / (mk + 2.0) * (vk_1 - vk) + 
#       max(1.0 - 100*sum((vk_1 - vk)^2) / sum((y[N1, 1:D] - y[N1 + 1, 1:D])^2), 
#           0.0) * rnorm(D) * 0.1
#     if(fn(y[n + 1, 1:D]) >= fn_old){
#       mk = mk + 1.0
#     }else{ 
#       mk = mk 
#     }
#     y[n + 1, D + 1] <- lp_f(y[n + 1, 1:D])
#     vk = vk_1;
#   }
#   y
# }

opt_path <- function(init, fn, gr, N1 = 30, N1_2 = 10, N_rep = 6, lp_f) {
  rho = 1.5; c = 2.0
  #init <- runif(D, -init_bound, init_bound)
  D <- length(init)
  y <- matrix(NA, nrow = N1 + N_rep*N1_2 + 1, ncol = D + 1)
  y[1, 1:D] <- init
  y[1, D + 1] <- lp_f(init)
  n0 = 1      # restart index
  for (j in 1:N_rep){
    cat("\n", j, "\n")
    if(j > 1.5){N = N1_2}else{N = N1}
    for (n in 1:N) {
      break_opt <- FALSE
      tryCatch(z <- optim(par = y[n0, 1:D],
                          fn = function(x) -fn(x),  # negate for maximization
                          gr = function(x) -gr(x),
                          method = "L-BFGS-B",
                          hessian = TRUE,
                          control = list(maxit = n, factr = 1e10#, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        print("Error in obtaining optimization path.")
        #return(y[1:(n0 + n), 1:(D + 1)])
        n = n - 1
        break
      }
      y[n0 + n, 1:D] <- z$par
      y[n0 + n, D + 1] <- lp_f(z$par)
      
      # break if no change in objective
      printf("n = %4d, last = %f;   this = %f",
             n0 + n , y[n + n0 - 1, D + 1], y[n + n0, D + 1])
      print(z$counts)
      if (y[n + n0 - 1, D + 1] == y[n + n0, D + 1]) {
        n = n - 1; break
      }
      if (abs(y[n + n0 - 1, D + 1] - y[n + n0, D + 1]) < 0.01 ) {
        break
      }
    }
    #next_opt <- FALSE
    # tryCatch(cholH <- chol(z$hessian),
    #          error = function(e) { 
    #            cat(": restart L-BFGS \n")
    #            next_opt <<- TRUE})
    # if(next_opt){ n0 = n0 + n; next}
    if (n == 0){
      diagH = diag(D)
      fn_last_max = fn(y[n0, 1:D])
    }else{
      diagH = abs(diag(z$hessian))
      fn_last_max = -z$value
    }
    n0 = n0 + n + 1
    u = matrix(rnorm(D*20), ncol = 20)
    #u2 = solve(cholH) %*% u %*% diag(1 / sqrt(colSums(u^2)))
    u2 = diag(1/sqrt(diagH))%*% u %*% diag(1 / sqrt(colSums(u^2)))
    
    fn_max = 0.0 
    max_ind = 0
    for(l in 1:ncol(u)){
      skip_flag = FALSE
      tryCatch(f_test <- fn(y[n0 - 1, 1:D] + c * rho^j * u2[, l]),
               error = function(e) { skip_flag <<- TRUE})
      if(skip_flag){next}
      if(max_ind == 0.0){fn_max = f_test; max_ind = l}
      else if(f_test > fn_max){
        fn_max = f_test; max_ind = l
      }
    }
    if((max_ind > 0) && (fn_max < (fn_last_max - 100*D))){
      # max_ind = 0;
      # y[n0, 1:D] = y[n0 - 1, 1:D]
      y[n0, 1:D] = y[n0 - 1, 1:D] + c * rho^j * u2[, max_ind]
    } else if (max_ind == 0){
      y[n0, 1:D] = y[n0 - 1, 1:D]
    } else{
      y[n0, 1:D] = y[n0 - 1, 1:D] + c * rho^j * u2[, max_ind]
    }
    y[n0, D + 1] = lp_f(y[n0, 1:D])
  }
  y[1:(n0-1), ]
}
# plot(30:(n0 - 1), y[30:(n0-1), D+1], type = "l")
# abline(h = c(INV))

# opt_path_2 <- function(x_init, fn, gr, N=40, lp_f) {
#   D <- length(x_init2)
#   y <- matrix(NA, nrow = N + 1, ncol = D + 1)
#   y[1, 1:D] <- x_init2
#   y[1, D + 1] <- fn(x_init2)
#   gt <- -gr(x_init2)
#   rho = 0.3
#   eps = 0.00001
#   Egt = sum(gt^2)
#   Edx = 0.0
#   for(n in 1:N){
#     gt <- -gr(y[n, 1:D])
#     Egt <- rho*Egt + (1 - rho)*sum(gt^2)
#     y[n + 1, 1:D] = y[n, 1:D] - sqrt(Edx + eps) / sqrt(Egt + eps) * gt #+ 
#       #max(0.5 - sqrt(Edx + eps) / sqrt(Egt + eps) * sqrt(sum(gt^2)), 0) * 
#       #(rnorm(D)*0.1)
#     cat(sqrt(Edx + eps) / sqrt(Egt + eps) * sqrt(sum(gt^2)), "\n")
#     y[n + 1, D + 1] = fn(y[n + 1, 1:D])
#     Edx = rho * Edx + (1 - rho)*sum((y[n + 1, 1:D] - y[n, 1:D])^2)
#   }
#   y
# }

opt_path_x <- function(init, fn, gr, N1 = 30, N_rep = 6, lp_f) {
  const = 2.0
  stepsize = 0.05
  #init <- runif(D, -init_bound, init_bound)
  D <- length(init)
  y <- matrix(NA, nrow = N1 * N_rep + 1, ncol = D + 1)
  y[1, 1:D] <- init
  y[1, D + 1] <- lp_f(init)
  n0 = 1      # restart index
  for (n in 1:N1) {
    break_opt <- FALSE
    tryCatch(z <- optim(par = y[1, 1:D],
                        fn = function(x) -fn(x),  # negate for maximization
                        gr = function(x) -gr(x),
                        method = "L-BFGS-B",
                        hessian = TRUE,
                        control = list(maxit = n, factr = 1e10#, ndeps = 1e-8 #,
                                       #trace = 6, REPORT = 1
                        )),
             error = function(e) { break_opt <<- TRUE})
    if(break_opt) {
      print("Error in obtaining optimization path.")
      #return(y[1:(n0 + n), 1:(D + 1)])
      n = n - 1
      break
    }

    y[n + 1, 1:D] <- z$par
    y[n + 1, D + 1] <- lp_f(z$par)

    # break if no change in objective
    printf("n = %4d, last = %f;   this = %f",
           n , y[n , D + 1], y[n + 1, D + 1])
    print(z$counts)
    if (y[n , D + 1] == y[n + 1, D + 1]) {
      n = n - 1; break
    }

    # if (abs(y[n + n0 - 1, D + 1] - y[n + n0, D + 1]) < 0.1 ) {
    #   break
    # }
  }
  # y[1:(n + 1), ]
  # Hessian <- z$hessian
  # next_opt <- FALSE
  # tryCatch(cholH <- chol(Hessian),
  #          error = function(e) { next_opt <<- TRUE})
  # if(next_opt){ n0 = n0 + n; next}
  n0 = n + 1
  for (j in 1:N_hmc){
    cat("\n", j, "\n")
    theta_up <- y[n0, 1:D]
    g0 <- gr(theta_up)
    lp0 <- fn(theta_up)
    rho0 =  const^j * (sqrt(length(theta_up)) / sqrt(sum(g0^2)))*g0
    rho <- rho0 + 0.5 * stepsize * g0
    break_flag <- FALSE
    draws <- c()
    for (l in 1:int_time){
      theta_up <- theta_up + stepsize * rho
      y[n0 + l, 1:D] <- theta_up
      break_flag <-  FALSE
      tryCatch(draws[l] <- fn(theta_up),
               error = function(e) { break_flag <<- TRUE })
      if(break_flag | is.na(draws[l])){l = l-1; break}
      y[n0 + l, D + 1] <- lp_f(theta_up)
      cat("\t", draws[l], "\t")
      if(l < int_time){
        tryCatch(g1 <- gr(theta_up),
                 error = function(e) { break })
        if(break_flag | any(is.na(g1))){break}
        rho <- rho + stepsize * g1
      }
    }
    l_max = which.max(draws)
    if(draws[l_max] > lp0){
      n0 = n0 + l_max
    }
  }
  y[1:n0, ]

  y[1:n0, D+1]

  n0 = n0 + n + 1
  u = matrix(rnorm(D*10), ncol = 10)
  u2 = solve(cholH) %*% u %*% diag(1 / sqrt(colSums(u^2)))
  fn_max = 0.0
  max_ind = 0.0
  for(l in 1:10){
    skip_flag = FALSE
    tryCatch(f_test <- fn(y[n0 - 1, 1:D] + rho^j * u2[, l]),
             error = function(e) { skip_flag <<- TRUE})
    if(skip_flag){next}
    if(max_ind == 0.0){fn_max = f_test; max_ind = l}
    else if(f_test > fn_max){
      fn_max = f_test; max_ind = l
    }
  }
  y[n0, 1:D] = y[n0 - 1, 1:D] + rho^j * u2[, max_ind]
  y[n0, D + 1] = lp_f(y[n0, 1:D])
  y[1:(n0-1), ]
}


opt_path_GDs <- function(x_init, fn, gr, N1 = 30, N1_2 = 10, N_rep = 6, lp_f) {
  rho = 1.1
  #x_init <- runif(D, -init_bound, init_bound)
  D <- length(x_init)
  y <- matrix(NA, nrow = N1 + 2, ncol = D + 1)
  y[1, 1:D] <- x_init
  y[1, D + 1] <- lp_f(x_init)
  n0 = 1      # restart index
  for (j in 1:N_rep){
    cat("\n", j, "\n")
    if(j > 1.5){N = N1_2}else{N = N1}
    for (n in 1:N) {
      break_opt <- FALSE
      tryCatch(z <- optim(par = y[n0, 1:D],
                          fn = function(x) -fn(x),  # negate for maximization
                          gr = function(x) -gr(x),
                          method = "L-BFGS-B",
                          hessian = TRUE,
                          control = list(maxit = n, factr = 1e10#, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        print("Error in obtaining optimization path.")
        #return(y[1:(n0 + n), 1:(D + 1)])
        n = n - 1
        break
      }
      
      y[n0 + n, 1:D] <- z$par
      y[n0 + n, D + 1] <- lp_f(z$par)
      
      # break if no change in objective
      printf("n = %4d, last = %f;   this = %f",
             n0 + n , y[n + n0 - 1, D + 1], y[n + n0, D + 1])
      print(z$counts)
      if (y[n + n0 - 1, D + 1] == y[n + n0, D + 1]) {
        n = n - 1; break
      }
      # if (abs(y[n + n0 - 1, D + 1] - y[n + n0, D + 1]) < 0.1 ) {
      #   break
      # }
    }
    #next_opt <- FALSE
    
    # tryCatch(cholH <- chol(z$hessian),
    #          error = function(e) { 
    #            cat(": restart L-BFGS \n")
    #            next_opt <<- TRUE})
    # if(next_opt){ n0 = n0 + n; next}
    diagH = abs(diag(z$hessian))
    n0 = n0 + n + 1
    u = matrix(rnorm(D*10), ncol = 10)
    #u2 = solve(cholH) %*% u %*% diag(1 / sqrt(colSums(u^2)))
    u2 = diag(1/sqrt(diagH))%*% u %*% diag(1 / sqrt(colSums(u^2)))
    fn_last_max = z$value
    fn_max = 0.0 
    max_ind = 0
    for(l in 1:10){
      skip_flag = FALSE
      tryCatch(f_test <- fn(y[n0 - 1, 1:D] + rho^j * u2[, l]),
               error = function(e) { skip_flag <<- TRUE})
      if(skip_flag){next}
      if(max_ind == 0.0){fn_max = f_test; max_ind = l}
      else if(f_test > fn_max){
        fn_max = f_test; max_ind = l
      }
    }
    if((max_ind > 0) && (fn_max < (fn_last_max - D))){
      max_ind = 0;
      y[n0, 1:D] = y[n0 - 1, 1:D]
    } else if (max_ind == 0){
      y[n0, 1:D] = y[n0 - 1, 1:D]
    } else{
      y[n0, 1:D] = y[n0 - 1, 1:D] + rho^j * u2[, max_ind]
    }
    y[n0, D + 1] = lp_f(y[n0, 1:D])
  }
  y[1:(n0-1), ]
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
opt_path_stan <- function(model, data, N1, N_rep, init_bound = 2) {
  # require chains > 0, iter > 0, so use Fixed_param to avoid work
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) log_prob(posterior, theta, adjust_transform = FALSE, 
                                 gradient = TRUE)[1]
  gr <- function(theta) grad_log_prob(posterior, theta, adjust_transform = FALSE)
  lp_f <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE, 
                                 gradient = TRUE)[1]
  out <- opt_path(init, fn, gr, N1, N1_2, N_rep, lp_f)
  return(out)
}
#opt_path_stan(model, data, N1, N_rep, init_bound = 2)

opt_path_stan_parallel <- function(seed_list, mc.cores, 
                                   model, data, N1, N1_2, N_rep, init_bound){
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
  out <- mclapply(init, opt_path, fn = fn, gr = gr, N1 = N1, N1_2 = N1_2,
                  N_rep = N_rep, lp_f = lp_f, mc.cores = mc.cores)
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
  draws <- c()
  rho <- rnorm(length(init_param_unc)) + 0.5 * stepsize * gr(theta_up)
  for (l in 1:int_time){
    theta_up <- theta_up + stepsize * rho
    draws[l] <- fn(theta_up)
    if(l < int_time){
      rho <- rho + stepsize * gr(theta_up)
    }
  }
  draws
}


increased <- function(lps) {
  # simple comparison of endpoint rather than whole path
  lps[length(lps)] > lps[1]
}

stuck <- function(lps) {
  # simple comparison of endpoint rather than whole path
  lps[length(lps)] == lps[1]
}

is_typical <- function(model, data, init_param_unc, M, int_time, lp_0) {
  increase_count <- 0
  stuck_count <- 0
  for (m in 1:M) {
    posterior <- to_posterior(model, data)
    init_fun <- function(chain_id) constrain_pars(posterior, init_param_unc)
    # find a stepsize #
    fit_0 <- sampling(model, data = data, init = init_fun,
                      chains = 1, iter = 2, warmup = 1, refresh = 0,
                      control = list(metric = "unit_e",
                                     adapt_engaged = TRUE,
                                     max_treedepth = 1),
                      save_warmup = TRUE)
    stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"]/2
    
    # generate hamiltonian dynamics
    fn <- function(theta) log_prob(posterior, theta, adjust_transform = TRUE, 
                                   gradient = TRUE)[1]
    gr <- function(theta) grad_log_prob(posterior, theta, 
                                        adjust_transform = FALSE)
    lps <- lp_draws(posterior, init_param_unc, int_time, fn, gr, stepsize)
    increase_count <- increase_count + sum(lps > lp_0) # increased(lps)
    #stuck_count <- stuck_count + (lps[length(lps)] == lp_0) #stuck(lps)
  }
  return(c(increase_count, stuck_count))
}

find_typical <- function(param_path, model, data, M = 4, int_time = 6) {
  typical_index <- c()          # return the index of sample in param_path that is identified as a good initial
  N <- dim(param_path)[1]
  D <- dim(param_path)[2] - 1   # includes objective in last position
  for (n in 1:N) {
    increase_counts <- is_typical(model, data, param_path[n, 1:D], M, int_time,
                                  param_path[n, D + 1])
    # printf("n = %3d;  increase proportion = %3.2f",
    #        n, increase_prop)
    printf("n = %3d;  increase count = %3d, stick count = %3d",
           n, increase_counts[1], increase_counts[2])
    #print(param_path[n, 1:(D + 1)], digits = 2)
    if(increase_counts[2] == M){next} # if all chains divergent, skip
    
    increase_prop = increase_counts[1] / (M * int_time - increase_counts[2])
    # declare typical if in central 90% interval of random increase/decrease
    # lb = qbinom(0.05, M, 0.5) / M
    # ub = qbinom(0.95, M, 0.5) / M
    lb = qbinom(0.1, (M * int_time), 0.5) / (M * int_time)
    ub = qbinom(0.9, (M * int_time), 0.5) / (M * int_time)
    if (increase_prop >= lb && increase_prop <= ub){
      print(param_path[n, 1:(D + 1)], digits = 2)
      typical_index <- c(typical_index, n)
    }
  }
  return(typical_index)
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

