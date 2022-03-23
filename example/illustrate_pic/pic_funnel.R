setwd("./posteriordb")
library(ggplot2)
library(posteriordb)
library(RColorBrewer)
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posteriordb)
library(posterior)
library(cmdstanr)
library(bayesplot)
library(gridExtra)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

## tuning parameters (default) ##
N1 = 1000    # maximum iters in optimization
factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
lmm = 6 # maximum histogram size

## check 2-d funnel ##
model <- stan_model(file = "../example/illustrate_pic/funnel.stan")
model_cmd <- cmdstan_model(stan_file = "../example/illustrate_pic/funnel.stan")

# set up lmm #
posterior <- to_posterior(model, data = NULL)
D <- get_num_upars(posterior)
cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")

set.seed(1)
init <- c(7, -8) #c(7, -7)
fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                gradient = TRUE)[1] 
gr <- function(theta) -grad_log_prob(posterior, theta, adjust_transform = TRUE)

# preallocation
sample_pkg <- list() # save the sample info for each mode
DIV_fit <- list() # save the divergence related info for each mode

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
    next}
  g1 <- gr(y[1, 1:D]) # record the gradient of initials
  if(any(is.na(g1))){
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
  if(LBFGS_fail){next}
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

# estimate DIV for all approximating Gaussians and save results of the one with maximum DIV
E <- rep(1, D)    # initial diag inv Hessian
Ykt_h <- NULL; Skt_h <- NULL # initialize matrics for storing history of updates
t <- proc.time()
for (iter in 1:Iter){
  cat(iter, "\t")
  inc_flag <- check_cond(Ykt[iter, ], Skt[iter, ])
  if(inc_flag){
    E <- Form_init_Diag(E, Ykt[iter, ], Skt[iter, ]) # initial estimate of diagonal inverse Hessian
    Ykt_h <- updateYS(Ykt_h = Ykt_h, Ykt = Ykt[iter, ], lmm) # update Y and S matrix
    Skt_h <- updateYS(Ykt_h = Skt_h, Ykt = Skt[iter, ], lmm)
  }
  ill_distr = FALSE
  tryCatch(
    # generate matrics for forming approximted inverse Hessian
    sample_pkg[[iter]] <- Form_N_apx(X[iter + 1, ], G[iter + 1, ], Ykt_h, Skt_h, E, lmm),
    error = function(e) { ill_distr <<- TRUE})
  if(ill_distr){ next }
  if(is.na(sample_pkg[1])){ next }
  DIV_fit[[iter]] <- est_DIV(sample_pkg[[iter]], N_sam = 100, fn, label = "ELBO") #lCDIV #lADIV  #lIKL #ELBO
} 
proc.time() -t 

# reference samples
fit_stan_flag <- FALSE
if(fit_stan_flag){
  fit_stan <- model_cmd$sample(
    chains=4, parallel_chains=4,
    seed = 1,
    adapt_delta = 0.9,
    thin = 100,
    iter_warmup = 50000,
    iter_sampling = 250000,
    show_messages = TRUE,
    sig_figs = 16)
  
  fit_stan$print("lp__")
  fit_stan$save_object(file = "../example/illustrate_pic/funnel_ref.RDS")
  p1 <- mcmc_trace(fit_stan$draws("lp__"), iter1 = 1) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  #' Check whether parameters have reasonable values
  draws <- fit_stan$draws()
  Sum_M <- summarise_draws(subset(draws, variable=c('log_sigma','alpha'), 
                                  regex=TRUE))
  summary(Sum_M$ess_bulk); summary(Sum_M$ess_tail)
  posterior <- to_posterior(model, data = NULL)
  
  draws <- as_draws_matrix(draws)
  Stan_draws <- matrix(draws, nrow = dim(draws)[1])
  colnames(Stan_draws) <- dimnames(draws)$variable
  unconstrained_draws <- unconstrain_cmd_draws(Stan_draws, posterior)
  colnames(unconstrained_draws) <- dimnames(draws)$variable[-1]
  
  lp_INV <- quantile(Stan_draws[, "lp__"], c(0.005, 0.995))
  
  save(file = "../example/illustrate_pic/reference_funnel.RData",
       list = c("unconstrained_draws", "lp_INV"))
}else{
  load("../example/illustrate_pic/reference_funnel.RData")
}

## plot figures ##
ref_dat <- data.frame(x = unconstrained_draws[, 2], y = unconstrained_draws[, 1])

alpha_g = seq(-20, 20, by = 0.1)
log_sigma_g = seq(-20, 20, by = 0.1)
locs <- expand.grid(log_sigma_g, alpha_g)
lp_ls <- -apply(locs, 1, fn)
lp_ls[lp_ls < -8] = NA

plot_list = list()
for(iter_trun in 3:Iter){
  cat(iter_trun, "\t")
  #iter_trun = 7 #4 5 6 7 8 9 
  density_lp_dat <- data.frame(x = locs[, 2], y = locs[, 1],
                               lp = lp_ls, n = as.numeric(!is.na(lp_ls)))
  opt_dat <- data.frame(optim_x = y[2:iter_trun, 2], 
                        optim_y = y[2:iter_trun, 1], 
                        optim_lp = y[2:iter_trun, 3])
  
  apx_dat <- data.frame(x = DIV_fit[[iter_trun-2]]$repeat_draws[2, ],
                        y = DIV_fit[[iter_trun-2]]$repeat_draws[1, ])
  
  p <- ggplot() +
    stat_density_2d(data = ref_dat, aes(x = x, y = y, fill = ..level..),
                    geom = "polygon", colour="#feb24c", size = 0.5) + 
    scale_fill_distiller(palette="YlOrBr", direction=1) +
    geom_contour(data = density_lp_dat, 
                 aes(x, y, z = lp, colour = after_stat(level)), na.rm = TRUE) + 
    scale_x_continuous(expand = c(0, 0), limits = c(-15, 15)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-13, 10)) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) + 
    guides(colour=guide_legend(title="log density")) +
    geom_point(data = opt_dat, 
               aes(x = optim_x, y = optim_y), size = 2) +
    geom_path(data = opt_dat, 
              aes(x = optim_x, y = optim_y)) +
    stat_ellipse(data = apx_dat, aes(x, y), size = 0.8) +
    xlab(paste0("iteration ", iter_trun - 1, "\n estimated ELBO: ", 
                round(DIV_fit[[iter_trun-2]]$DIV, digits = 1))) + 
    ylab("") + theme_bw() + theme(legend.position='none')
  
  print(p)
  plot_list[[iter_trun - 2]] <- p
  ggsave(paste0("funnel", iter_trun - 1, ".eps"), #"8-school_opt_tr22.eps"
         plot = p,
         device = cairo_ps,
         path = "../example/illustrate_pic/pic",
         width = 5.0, height = 5.0, units = "in")
  
}

p_comb <- grid.arrange(plot_list[[2]], plot_list[[3]], plot_list[[4]], 
                       plot_list[[5]], plot_list[[8]], plot_list[[12]],
                       ncol = 3, nrow = 2)

ggsave(paste0("funnel_example.png"), #"8-school_opt_tr22.eps"
       plot = p_comb,
       path = "../example/illustrate_pic/pic",
       width = 9.0, height = 6.0, units = "in")


ggsave(paste0("funnel_example.eps"), #"8-school_opt_tr22.eps"
       plot = p_comb,
       path = "../example/illustrate_pic/pic",
       width = 9.0, height = 6.0, units = "in")

