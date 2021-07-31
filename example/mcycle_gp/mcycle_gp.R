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
model <- stan_model(file = "../example/mcycle_gp/mcycle_gp_accel_gp.stan")

### load data ###
pd <- pdb_local()
po <- posterior("mcycle_gp-accel_gp", pdb = pd)
data <- get_data(po)

fit_pf_flag <- FALSE 
# preallocate results #
if(fit_pf_flag){
  mc.cores = parallel::detectCores() - 2
  ## tuning parameters
  init_bound = 2.0 # parameter for initial distribution 
  N1 = 1000    # maximum iters in optimization
  factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
  N_sam_DIV = 5   # samples for ELBO evaluation
  N_sam = 100
  lmm = 100 # history size
  seed_list = 1:20
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  t <- proc.time()
  opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data,
                                  init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                  factr_tol, lmm) # plot for 8school init_bound = 15
  print(proc.time() - t)
  
  pick_samples_IR <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)
  #pick_samples_IR <- Imp_Resam_WOR(opath, n_inits = 100, seed = 1)
  #pick_samples_IR <- filter_samples_resam(opath, n_inits = 100, seed = 1)
  pick_samples <- pick_samples_IR
  pf_fn_calls <- sapply(opath, f <- function(x){x$fn_call})
  pf_gr_calls <- sapply(opath, f <- function(x){x$gr_call})
  save(file = "../example/mcycle_gp/opath_org100h.RData",
       list = c("opath", "pick_samples", "pf_fn_calls", "pf_gr_calls"))
  
}else{
  load("../example/mcycle_gp/opath_org100h.RData")
}

## load reference samples ##
gsd <- reference_posterior_draws(po)
posterior <- to_posterior(model, data)
unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                    unconstrained_draws[[3]], unconstrained_draws[[4]],
                    unconstrained_draws[[5]], unconstrained_draws[[6]],
                    unconstrained_draws[[7]], unconstrained_draws[[8]],
                    unconstrained_draws[[9]], unconstrained_draws[[10]])

## check plots ##
check_dim <- c(1, 2)
for(first_check in 1:20){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  
  unconstrained_draws <- ref_samples 
  
  plot(unconstrained_draws[, check_dim[1]], unconstrained_draws[, check_dim[2]], 
       col = "grey", 
       xlim = range(unconstrained_draws[, check_dim[1]],
                    pick_samples[check_dim[1], ]),
       ylim = range(unconstrained_draws[, check_dim[2]],
                    pick_samples[check_dim[2], ]),
       main = paste0("element ", check_dim[1], " and ", check_dim[2]))
  
  points(x = pick_samples[check_dim[1], ], y = pick_samples[check_dim[2], ], 
         col = "orange", pch = 16)
  
  # text(x = pick_samples[check_dim[1], ], y = pick_samples[check_dim[2], ], 
  #      labels=round(lps_samples, digits = 1), cex=0.9, font=2, col = "blue")
  
  readline(prompt="Press [enter] to continue:")
}

#(1, 2) (105, 106) (115, 116)
check_dim <- c(45, 46)
pf_lp <- apply(pick_samples, 2, f <- function(x){
  log_prob(posterior, x, adjust_transform = TRUE, gradient = TRUE)[1]
})
table(pf_lp)
#summary(Stan_draws[, "lp__"])
#quantile(Stan_draws[, "lp__"], c(0.05, 0.95))


## check the Wasserstein distance ##
library(transport)
ref_samples = unconstrained_draws

# pathfinder #
if(ncol(pick_samples) == 1){
  a = wpp(rbind(t(pick_samples), t(pick_samples)), 
          mass = rep(1 / 2, 2))
}else{
  a = wpp(t(pick_samples), 
          mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
}
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_pf <- wasserstein(a, b, p = 2); w_d_pf #19.06423 14.99645 #19.849  #19.36716 h100


## cost ##
# pf
sapply(opath, f <- function(x){x$fn_call})
sapply(opath, f <- function(x){x$gr_call})


## optimization ##
# opt <- mod$optimize(data=data, init=2, seed = 2) #2 321
# Est_mu_opt<- summarise_draws((subset(opt$draws(), variable = "mu", regex=TRUE)), "mean")
# Est_sigma_opt<- summarise_draws((subset(opt$draws(), variable = "sigma", regex=TRUE)), "mean")
# tmp <- rstan:::create_skeleton(posterior@model_pars, posterior@par_dims)
# skeleton <- tmp[names(tmp) != "lp__"]
# constrained_s <- rstan:::rstan_relist(opt$draws()[-1], skeleton)
# unconstrained_opt <-posterior@.MISC$stan_fit_instance$unconstrain_pars(constrained_s)


## check expectation ##
# expectation based on reference samples #
constrained_ref_sams <- 
  apply(ref_samples, 1, f <- function(x){constrain_pars(posterior, x)$mu})
Est_mu <- rowMeans(constrained_ref_sams)

# expectation based on Pathfinder approximate draws #
constrained_sams <- 
  apply(pick_samples, 2, f <- function(x){constrain_pars(posterior, x)$mu})
Est_mu_pf <- rowMeans(constrained_sams)

plot(Est_mu, Est_mu_pf, xlab = "E mu from reference samples",
     ylab = "E mu from Pathfinder samples")
abline(a = 0, b = 1)

dat_compar <- data.frame(Est_mu = Est_mu, Est_mu_pf = Est_mu_pf)
p_compar <- ggplot(dat_compar, aes(x=Est_mu, y=Est_mu_pf)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype="dashed") + 
  theme_bw() +
  labs(x="Expected mean (reference)", y="Expected mean (Pathfinder)")
p_compar
ggsave("mcycle_pred_compar.eps", #"8-school_opt_tr22.eps"
       plot = p_compar,
       device = cairo_ps,
       path = "../example/mcycle_gp/",
       width = 4.0, height = 4.0, units = "in")


plot(Est_mu, Est_mu_opt$mean, xlab = "E mu from reference samples",
     ylab = "E mu from Pathfinder samples")
abline(a = 0, b = 1)


# expectation based on reference samples #
constrained_ref_sams <- 
  apply(ref_samples, 1, f <- function(x){constrain_pars(posterior, x)$sigma})
Est_sigma <- rowMeans(constrained_ref_sams)

# expectation based on Pathfinder approximate draws #
constrained_sams <- 
  apply(pick_samples, 2, f <- function(x){constrain_pars(posterior, x)$sigma})
Est_sigma_pf <- rowMeans(constrained_sams)


plot(Est_sigma, Est_sigma_pf, 
     xlab = "E sigma from reference samples",
     ylab = "E sigma from Pathfinder samples")
abline(a = 0, b = 1)


# ## for fun ##
# init_ls <- list()
# for (s in 1:length(seed_list)){
#   init_ls[[s]] <- ref_samples[s,]
# }
# t <- proc.time()
# opath <- opt_path_stan_init_parallel(
#   init_ls, mc.cores, model, data = data, init_bound = 2.0, 
#   N1, N_sam_DIV, N_sam, factr_tol, lmm, seed_list)
# print(proc.time() - t)
# pick_samples <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)


## plots ##
opt_tr_res <- get_opt_tr(opath)
check_dim <- c(1, 2)

dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])

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
  scale_x_continuous(expand = c(0, 0))+#, limits = c(-13, 13)) +
  scale_y_continuous(expand = c(0, 0))+#, limits = c(-40, 13)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
ggsave("mcycle_1_2_points.eps", 
       plot = p_check,
       device = cairo_ps,
       path = "../example/mcycle_gp/",
       width = 5.0, height = 5.0, units = "in")

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0))+#, limits = c(-13, 13)) +
  scale_y_continuous(expand = c(0, 0))+#, limits = c(-40, 13)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.5), size = 1) +
  geom_path(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.5)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("mcycle_1_2_tr.eps", #"8-school_opt_tr22.eps"
       plot = p_check1,
       device = cairo_ps,
       path = "../example/mcycle_gp/",
       width = 5.0, height = 5.0, units = "in")

# amplified plots
p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-4, 2))+#, limits = c(-4, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 6))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_path(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
             size = 3)+
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
print(p_check1)
ggsave("32_opt_tr_12_s_h100.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")




