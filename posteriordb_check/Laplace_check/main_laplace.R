rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
#options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
library(loo)
library(transport)
library(numDeriv)   # Hessian approximation
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local(path = "./") # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# Laplace tuning parameters #
init_bound = 2.0 # parameter for initial distribution 
N1 = 1000    # maximum iters in optimization
factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
N_sam = 100
lmm = 6 # maximum histogram size

load(file = "../results/lp_posteriordb_LBFGS_h6.RData")

pick_ind <- c(1:7, 13, 15:21, 31:32, 46, 47:48) # index to select the models in the summary
model_record <- c(1, 2, 3, 4, 6, 7, 8, 14, 20, 21, 23, 24, 25, 26, 27, 40,
                  41, 58, 61, 94)

# precompile all models
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
}

load("../results/wasserstein_100_default_W1_WR_updat.RData") # load wasserstein distance
Laplace_draw <- list()
Laplace_center <- list()
Laplace_cov <- list()
w_d_laplace <- c()
calls_lp_laplace <- array(data = 0, dim = c(length(model_record)))
calls_gr_laplace <- array(data = 0, dim = c(length(model_record)))

ELBO_Laplace = array(data = NA, dim = c(1, length(model_record)))
iter_Laplace = array(data = NA, dim = c(1, length(model_record)))
opt_fit_ls <- list()
check_maxit_hit = rep(0, length(model_record))
## check laplace ##
#load("../results/Laplace_results.RData")
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  # set up initials
  data = get_data(po)
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1] 
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  
  seed = 0
  LBFGS_fail <- TRUE
  N_try <- 0
  while(LBFGS_fail){   # if L-BFGS cannot run with the initial, find another inital
    seed = seed + 1
    LBFGS_fail <- FALSE
    if(N_try == 0){
      N_try = N_try + 1
      set.seed(seed)
      init <- runif(D, -init_bound, init_bound)
    }else{
      print("\n reinitialize \n")
      N_try = N_try + 1
      set.seed(seed)
      init <- runif(D, -init_bound, init_bound)
    }
    tryCatch(y_1 <- -fn(init), 
             error = function(e) { LBFGS_fail <<- TRUE})
    if(LBFGS_fail | is.infinite(y_1)){
      LBFGS_fail <- TRUE
      calls_lp_laplace[i] = calls_lp_laplace[i] + 1  # record the evaluation of log-density of ill initials
      next}
    g1 <- gr(init) # record the gradient of initials
    if(any(is.na(g1))){
      calls_gr_laplace[i] = calls_gr_laplace[i] + 1  # record the evaluation of gradient of ill initials
      LBFGS_fail <- TRUE
      next
    }
    tryCatch(
      my_data <- capture.output(
        opt_fit_ls[[i]] <- optimx(par = init,
                     fn = fn,  # negate for maximization
                     gr = gr,
                     method = "L-BFGS-B",
                     control = list(maxit = ifelse((i == 5) | (i == 32), 
                                                   2000, N1), 
                                    pgtol = 0.0, 
                                    factr = factr_tol,
                                    trace = 6, REPORT = 1, lmm = lmm)),
        type = "output"), 
      error = function(e) { LBFGS_fail <<- TRUE})
    if(LBFGS_fail){ # fail the section of code that does the checking
      calls_lp_laplace[i] = calls_lp_laplace[i] + 1
      calls_gr_laplace[i] = calls_gr_laplace[i] + 1
      next}
    if(opt_fit_ls[[i]]$convcode == 9999){
      LBFGS_fail <- TRUE
      splited_output = unlist(lapply(my_data, f <- function(x){
        strsplit(as.character(x),split = " ")}))
      
      line_search_ind = which(splited_output == "SEARCH") + 1 # check line seach times the 
      
      eval_count <- sum(as.numeric(splited_output[line_search_ind])) + 
        length(line_search_ind) + 1
      calls_lp_laplace[i] = calls_lp_laplace[i] + eval_count  
      calls_gr_laplace[i] = calls_gr_laplace[i] + eval_count
    }
    if(opt_fit_ls[[i]]$convcode == 1){
      check_maxit_hit[i] <- 1
    }
  }
  cat("Hessian", "\t")
  
  Hess_est <- jacobian(gr, unlist(opt_fit_ls[[i]][1:D]))
  #Hess_est <- hessian(fn, unlist(opt_fit_ls[[i]][1:D]), method="Richardson")
  Laplace_center[[i]] <- unlist(opt_fit_ls[[i]][1:D])
  post_Hess <- TRUE
  tryCatch(
    Chol_Hess_est <- chol(Hess_est),
    error = function(e) { post_Hess <<- FALSE}
  )
  if(!post_Hess){
    cat("Hessian not positive definited")
    next
  }
  Laplace_cov[[i]] <- chol2inv(Chol_Hess_est)
  cat("sample", "\t")
  U = backsolve(Chol_Hess_est, 
                Matrix(rnorm(D * N_sam), nrow = D, ncol = N_sam)) + 
    unlist(opt_fit_ls[[i]][1:D])
  
  Laplace_draw[[i]] <- t(U)
  
}

# save(file = "../results/Laplace_results.RData",
#      list = c("Laplace_draw"))




## compute wasserstein distance ##
dimension_record <- c()
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  # set up initials
  data = get_data(po)
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  dimension_record[i] <- D
  
  if(is.null(Laplace_draw[[i]])){next}
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]],
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  # last samples of phase I #
  a_laplace = wpp(Laplace_draw[[i]],
                  mass = rep(1 / nrow(Laplace_draw[[i]]),
                             nrow(Laplace_draw[[i]])))
  b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
  w_d_laplace[i] <- wasserstein(a_laplace, b, p = 1) 
  print(summary(W_d_100_pf[, pick_ind[i]]))
  cat(w_d_laplace[i], "\n")
}

save(file = "../results/Laplace_results.RData",
     list = c("Laplace_draw", "w_d_laplace", "Laplace_center", 
              "Laplace_cov", "dimension_record", "pick_ind"))


# Code for Checking #
for(i in 1:length(model_record)){
  if(is.null(Laplace_draw[[i]])){next}
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  # set up initials
  data = get_data(po)
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  
  cat("benchmark lp INV:", round(lp_INV[, i], digits = 3), 
      "benchmark lp_ mean:", round(lp_mean[i], digits = 3), "\n")
  cat("dimension: ", D, "\n")
  print(summary(W_d_100_ADVI_mf[, pick_ind[i]]))
  print(summary(W_d_100_ADVI_fr[, pick_ind[i]]))
  print(summary(W_d_100_pf[, pick_ind[i]]))
  print(summary(W_d_100_pf_IR[, pick_ind[i]]))
  cat(w_d_laplace[i], "\n")
}

load("../results/multi_pf_samples_default.RData")

dim_check = c(5, 3)
plot(ref_samples[, dim_check[1]], ref_samples[, dim_check[2]], 
     main = modelname, xlab = paste("element", dim_check[1]), 
     ylab = paste("element", dim_check[2]))
# points(x = Laplace_center[[i]][dim_check[1]], 
#        y = Laplace_center[[i]][dim_check[2]], col = "orange")
# points(x = Laplace_draw[[i]][, dim_check[1]], 
#        y = Laplace_draw[[i]][, dim_check[2]], col = "orange")
# points(x = lp_opath[[i]]$opath[[12]]$DIV_save$repeat_draws[dim_check[1], ],
#        y = lp_opath[[i]]$opath[[12]]$DIV_save$repeat_draws[dim_check[2], ],
#        col = "yellow")
points(x = lp_multi_opath[[i]][[1]]$opath[[3]]$y[, dim_check[1]],
       y = lp_multi_opath[[i]][[1]]$opath[[3]]$y[, dim_check[2]],
  col = "yellow")
points(x = lp_multi_opath[[i]][[1]]$pick_samples[dim_check[1], ],
       y = lp_multi_opath[[i]][[1]]$pick_samples[dim_check[2], ], col = "green")




