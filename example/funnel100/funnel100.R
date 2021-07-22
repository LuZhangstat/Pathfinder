rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(bayesplot)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posteriordb)
library(posterior)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)

### generate reference posterior samples ###

# compile the model
model <- stan_model(file = "../example/funnel100/funnel100.stan")
model_cmd <- cmdstan_model(stan_file = "../example/funnel100/funnel100.stan")

### load data ###
## generate reference posterior samples ##
fit_stan_flag <- FALSE
if(fit_stan_flag){
  fit_stan <- model_cmd$sample(
                        chains=4, parallel_chains=4,
                        seed = 1,
                        adapt_delta = 0.9,
                        thin = 300,
                        iter_warmup = 100000,
                        iter_sampling = 750000,
                        show_messages = TRUE,
                        sig_figs = 16)
  
  fit_stan$print("lp__")
  fit_stan$save_object(file = "../example/funnel100/funnel100_ref.RDS")
  p1 <- mcmc_trace(fit_stan$draws("lp__"), iter1 = 1) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  #' Check whether parameters have reasonable values
  draws <- fit_stan$draws()
  Sum_M <- summarise_draws(subset(draws, variable=c('log_sigma','alpha'), 
                         regex=TRUE))
  summary(Sum_M$ess_bulk); summary(Sum_M$ess_tail)

}

## save the unconstrained reference posterior samples ##
posterior <- to_posterior(model, data = NULL)

## load reference posterior samples ##
if(fit_stan_flag){
  draws <- as_draws_matrix(draws)
  Stan_draws <- matrix(draws, nrow = dim(draws)[1])
  colnames(Stan_draws) <- dimnames(draws)$variable
  unconstrained_draws <- unconstrain_cmd_draws(Stan_draws, posterior)
  colnames(unconstrained_draws) <- dimnames(draws)$variable[-1]
  
  lp_INV <- quantile(Stan_draws[, "lp__"], c(0.005, 0.995))
  
  save(file = "../example/funnel100/reference_funnel100.RData",
       list = c("unconstrained_draws", "lp_INV"))
}else{
  load("../example/funnel100/reference_funnel100.RData")
}

## Pathfinder ##
fit_pf_flag <- FALSE
if(fit_pf_flag){
  mc.cores = parallel::detectCores() - 2
  ## tuning parameters
  init_bound = 2.0 # parameter for initial distribution 
  N1 = 1000    # maximum iters in optimization
  factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
  N_sam_DIV = 5   # samples for ELBO evaluation
  N_sam = 100
  lmm = 6 # history size
  seed_list = 1:20
  
  D <- get_num_upars(posterior)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  t <- proc.time()
  opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data = NULL,
                                  init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                  factr_tol, lmm) # plot for 8school init_bound = 15
  print(proc.time() - t)
  
  pick_samples_IR <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)
  pick_samples <- pick_samples_IR
  pf_fn_calls <- sapply(opath, f <- function(x){x$fn_call})
  pf_gr_calls <- sapply(opath, f <- function(x){x$gr_call})
  save(file = "../example/funnel100/opath.RData",
       list = c("opath", "pick_samples", "pf_fn_calls", "pf_gr_calls"))
}else{
  load("../example/funnel100/opath.RData")
}



## check plots ##
check_dim <- c(1, 2)
for(first_check in 6:10){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  
  plot(unconstrained_draws[, check_dim[1]], unconstrained_draws[, check_dim[2]], 
       col = "grey", 
       xlim = c(-10, 10),
         # range(unconstrained_draws[, check_dim[1]],
         #            pick_samples[check_dim[1], ]),
       ylim = c(-10, 10),
         # range(unconstrained_draws[, check_dim[2]],
         #            pick_samples[check_dim[2], ]),
       main = paste0("element ", check_dim[1], " and ", check_dim[2]),
       xlab = colnames(unconstrained_draws)[check_dim[1]], 
       ylab = colnames(unconstrained_draws)[check_dim[2]])
  
  points(x = pick_samples[check_dim[1], ], y = pick_samples[check_dim[2], ], 
         col = "orange", pch = 16)
  
  readline(prompt="Press [enter] to continue:")
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


opt_tr_res <- get_opt_tr(opath)#opath #lp_opath[[15]]$opath
check_dim <- c(2, 1)

dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = unconstrained_draws[, check_dim[1]],
                        y = unconstrained_draws[, check_dim[2]])

dta_opt <- data.frame(
  optim_x = opt_tr_res$opt_tr[, check_dim[1]],
  optim_y = opt_tr_res$opt_tr[, check_dim[2]],
  optim_ind = opt_tr_res$ind_tr,
  tr_id = opt_tr_res$tr_id
)
dta_opt$tr_id = factor(dta_opt$tr_id)

p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-13, 13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-40, 13)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
ggsave("funnel100_points.eps", #"8-school_opt_tr22.eps"
       plot = p_check,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-13, 13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-40, 13)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.5), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.5)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("funnel100_tr.eps", #"8-school_opt_tr22.eps"
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")



plot(unconstrained_draws[, check_dim[1]], unconstrained_draws[, check_dim[2]], 
     col = "grey", 
     xlim = range(unconstrained_draws[, check_dim[1]],
                  pick_samples[check_dim[1], ]),
     ylim = range(unconstrained_draws[, check_dim[2]],
                  pick_samples[check_dim[2], ]),
     main = paste0("element ", check_dim[1], " and ", check_dim[2]),
     xlab = colnames(unconstrained_draws)[check_dim[1]], 
     ylab = colnames(unconstrained_draws)[check_dim[2]])

points(x = pick_samples[check_dim[1], ], y = pick_samples[check_dim[2], ], 
       col = "orange", pch = 16)

## for fun ##
init_ls <- list()
for (s in 1:length(seed_list)){
  init_ls[[s]] <- opath[[s]]$y[1, 1:D] + 10
}
t <- proc.time()
opath <- opt_path_stan_init_parallel(
  init_ls, mc.cores, model, data = NULL, init_bound = 2.0, 
  N1, N_sam_DIV, N_sam, factr_tol, lmm, seed_list)
print(proc.time() - t)
opath <- opath[-c(7, 13, 14)]
pick_samples <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)


