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

pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)

modelname <- pn[62]
printf("model: %s", modelname)
# pick model
po <- posterior(modelname, pdb = pd)

# compile the model
sc <- stan_code(po)
model <- stan_model(model_code = sc)
data = get_data(po)
posterior <- to_posterior(model, data)
write_stan_file(sc, dir = paste0(getwd(), "/modelcode"),
                basename = paste0(modelname, ".stan"))
file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
model_cmd <- cmdstan_model(file)

### load data ###
## generate reference posterior samples ##
fit_stan_flag <- FALSE
if(fit_stan_flag){
  fit_stan <- model_cmd$sample(data = data,
                               chains=4, parallel_chains=4,
                               seed = 1,
                               adapt_delta = 0.9,
                               thin = 10,
                               iter_warmup = 10000,
                               iter_sampling = 25000,
                               show_messages = TRUE,
                               sig_figs = 16)
  
  fit_stan$print("lp__")
  fit_stan$save_object(file = "../example/overian/overian_ref.RDS")
  p1 <- mcmc_trace(fit_stan$draws("lp__"), iter1 = 1) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  #' Check whether parameters have reasonable values
  draws <- fit_stan$draws()
  Sum_M <- summarise_draws(subset(draws, regex=TRUE))
  summary(Sum_M$ess_bulk); summary(Sum_M$ess_tail)
  
  draws <- as_draws_matrix(draws)
  Stan_draws <- matrix(draws, nrow = dim(draws)[1])
  colnames(Stan_draws) <- dimnames(draws)$variable
  unconstrained_draws <- unconstrain_cmd_draws(Stan_draws, posterior)
  
  lp_INV <- quantile(Stan_draws[, "lp__"], c(0.005, 0.995))
  
  save(file = "../example/overian/reference_overian.RData",
       list = c("unconstrained_draws", "lp_INV"))
}else{
  load("../example/overian/reference_overian.RData")
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
  lmm = 6     # history size
  seed_list = 1:20
  
  D <- get_num_upars(posterior)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  t <- proc.time()
  opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data,
                                  init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                  factr_tol, lmm) # plot for 8school init_bound = 15
  print(proc.time() - t)
  
  #pick_samples_IR <- filter_samples_resam(opath, n_inits = 1000, seed = 1)
  pick_samples <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)
  pf_fn_calls <- sapply(opath, f <- function(x){x$fn_call})
  pf_gr_calls <- sapply(opath, f <- function(x){x$gr_call})
  save(file = "../example/overian/opath.RData",
       list = c("opath", "pick_samples", "pf_fn_calls", "pf_gr_calls"))
}else{
  load("../example/overian/opath.RData")
}

## pathfinder with initials from reference samples##
fit_pf_init_flag = FALSE
if(fit_pf_init_flag){
  mc.cores = parallel::detectCores() - 2
  ## tuning parameters
  N1 = 1000    # maximum iters in optimization
  factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
  N_sam_DIV = 5   # samples for ELBO evaluation
  N_sam = 100
  lmm = 6 # history size
  seed_list = 1:20
  
  D <- get_num_upars(posterior)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  init_ls <- list()
  for (s in 1:length(seed_list)){
    init_ls[[s]] <- unconstrained_draws[s,]
  }
  t <- proc.time()
  opath <- opt_path_stan_init_parallel(
    init_ls, mc.cores, model, data = data, init_bound = 2.0, 
    N1, N_sam_DIV, N_sam, factr_tol, lmm, seed_list)
  print(proc.time() - t)
  pick_samples <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)
  pf_fn_calls <- sapply(opath, f <- function(x){x$fn_call})
  pf_gr_calls <- sapply(opath, f <- function(x){x$gr_call})
  save(file = "../example/overian/opath_init.RData",
       list = c("opath", "pick_samples", "pf_fn_calls", "pf_gr_calls"))
}

## check plots ##
check_dim <- c(1, 2)
check_dim <- c(1537, 1538)
check_dim <- c(1, 1538)
check_dim <- c(3074, 3075)
check_dim <- c(1631, 2655) #93, 1117, 1491 
check_dim <- c(1631, 3029) 
check_dim <- c(1, 1492)
check_dim <- c(530, 2068)
for(first_check in 1:10){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  
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
  
  readline(prompt="Press [enter] to continue:")
}


opt_tr_res <- get_opt_tr(opath) 
check_dim <- c(1, 3029)#c(1631, 3029)#c(3029, 2655) #c(1, 2)
#check_dim <- c(1, 1538)

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

# p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#   scale_fill_distiller(palette=4, direction=-1) +
#   scale_x_continuous(expand = c(0, 0), limits = c(-10, 20)) + #c(-40, 50)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(-30, 5)) + #c(-10, 18)) +
#   #geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
#   xlab("") + ylab("") +
#   theme(
#     legend.position='none'
#   )
# 
# p_check
# ggsave("overian_1_2.eps", #"8-school_opt_tr22.eps"
#        plot = p_check,
#        device = cairo_ps,
#        path = "../example/overian/",
#        width = 5.0, height = 5.0, units = "in")

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits =c(-40, 50)) + #c(-10, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-10, 18)) + #c(-30, 5)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.5), size = 1) +
  geom_path(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.5)) +
  scale_color_gradient(low="white", high="orange") +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("overian_1_2_tr.eps", #_init 
       plot = p_check1,
       device = cairo_ps,
       path = "../example/overian/",
       width = 5.0, height = 5.0, units = "in")


check_dim <- c(1, 2)
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

## check the prediction ##
# expectation based on reference samples #
fit_stan <- readRDS("../example/overian/overian_ref.RDS")
draws <- fit_stan$draws()
Sum_lh<- 
  summarise_draws(exp(subset(draws, variable = "log_lik", regex=TRUE)), "mean")
Sum_f<- 
  summarize_draws(plogis(subset(draws, variable = "f", regex=TRUE)), "mean")


# expectation based on Pathfinder approximate draws #
constrained_sams <- 
  apply(pick_samples, 2, f <- function(x){constrain_pars(posterior, x)$log_lik})
Est_lh <- rowMeans(exp(constrained_sams))
f_sams_pf <- 
  apply(pick_samples, 2, f <- function(x){constrain_pars(posterior, x)$f})
Est_f <- rowMeans(plogis(f_sams_pf))


plot(Sum_lh$mean, Est_lh, xlim = c(0, 1), ylim = c(0, 1), xlab = "ref sample",
     ylab = "1 ref sample", 
     main = "expected probability")
abline(a = 0, b = 1)

plot(Sum_f$mean, Est_f, xlab = "HMC", xlim = c(0, 1), ylim = c(0, 1),
     ylab = "Pathfinder")
abline(a = 0, b = 1)

dat_compar <- data.frame(Est_p = Sum_f$mean, Est_p_pf = Est_f)
p_compar <- ggplot(dat_compar, aes(x=Est_p, y=Est_p_pf)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype="dashed") + 
  theme_bw() + xlim(0, 1) + ylim(0, 1) + 
  labs(x="Probability (reference)", y="Probability (Pathfinder)")
p_compar
ggsave("overian_pred_compar.eps",  #_init
       plot = p_compar,
       device = cairo_ps,
       path = "../example/overian/",
       width = 4.0, height = 4.0, units = "in")



# check for fun #
# which.max(sapply(opath, f <-
#                    function(x){mean(-x$DIV_save$fn_draws - x$DIV_save$lp_approx_draws)}))
# pick_samples <- opath[[10]]$DIV_save$repeat_draws
# library(ggcorrplot)
# cor_M <- cor(unconstrained_draws)
# p_cor <- ggcorrplot(cor_M)

# library(transport)
# if(ncol(pick_samples) == 1){
#   a = wpp(rbind(t(pick_samples), t(pick_samples)),
#           mass = rep(1 / 2, 2))
# }else{
#   a = wpp(t(pick_samples),
#           mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
# }
# b = wpp(unconstrained_draws, mass = rep(1 / nrow(unconstrained_draws), 
#                                         nrow(unconstrained_draws)))
# w_d_pf <- wasserstein(a, b, p = 1)#wasserstein(a, b, p = 2); w_d_pf
# w_d_pf   #96.19 vs #98.15165


library("scatterplot3d")
#setEPS()
png(filename = "../example/overian/ovarian_mode.png", width = 600, 
    height = 400)
stp3 <- scatterplot3d(unconstrained_draws[, c(1631, 2655, 3029)], 
                      cex.symbols = 0.05, xlab = "", #bquote(lambda[93]), 
                      ylab = "", #bquote(lambda[1117]), 
                      zlab = "", #bquote(lambda[1491])
                      xlim = c(-8, 10),
                      tick.marks=FALSE, label.tick.marks=FALSE, angle = 45,
                      color = alpha("black", 0.4),
                      grid=TRUE, box=FALSE)

dev.off()

## check for fun ##
stp3 <- scatterplot3d(unconstrained_draws[, c(1631, 3021, 3029)], 
                      cex.symbols = 0.1, xlab = bquote(lambda[93]),
                      ylab = bquote(lambda[1483]), zlab = bquote(lambda[1491]))


stp3$points3d(t(pick_samples[c(1631, 2655, 3029), ]), col = "red", pch = 16)
check_dims <- c(1, 1631, 1538) #c(1631, 3021, 3029) #
par(mfrow=c(2,2))
for(ind in 1:20){
  est_ELBO = round(mean(-opath[[ind]]$DIV_save$fn_draws -
                          opath[[ind]]$DIV_save$lp_approx_draws), digits = 2)
  cat(ind, "\t", opath[[ind]]$DIV_save$DIV, "\t",
      mean(-opath[[ind]]$DIV_save$fn_draws -
             opath[[ind]]$DIV_save$lp_approx_draws), "\n") #c(1, 3029, 1631) #c(1631, 2655, 3029)
  pick_samples <- opath[[ind]]$DIV_save$repeat_draws
  stp3 <- scatterplot3d(unconstrained_draws[, check_dims], 
                        cex.symbols = 0.1,
                        xlim = range(c(opath[[ind]]$y[, check_dims[1]], 
                                       unconstrained_draws[, check_dims[1]])),
                        ylim = range(c(opath[[ind]]$y[, check_dims[2]], 
                                       unconstrained_draws[, check_dims[2]])),
                        zlim = range(c(opath[[ind]]$y[, check_dims[3]], 
                                       unconstrained_draws[, check_dims[3]])),
                        xlab = "", ylab = "", zlab = "",
                        angle = 45,
                        main = paste0(ind, "th pf, est ELBO:", est_ELBO))
  stp3$points3d(t(pick_samples[check_dims, ]), col = "red", pch = 16)
  stp3$points3d(opath[[ind]]$y[1:as.integer(nrow(opath[[ind]]$y) / 2), 
                               check_dims], col = "orange", pch = 16)
  stp3$points3d(opath[[ind]]$y[as.integer(nrow(opath[[ind]]$y) / 2):nrow(opath[[ind]]$y), 
                               check_dims], col = "green", pch = 16)
  stp3$points3d(t(opath[[ind]]$y[1, check_dims]), col = "yellow", pch = 16)
  stp3$points3d(t(opath[[ind]]$y[nrow(opath[[ind]]$y), check_dims]),
                col = "blue", pch = 16)
  readline(prompt="Press [enter] to continue:")
  f_sams_pf <- 
    apply(pick_samples, 2, f <- function(x){constrain_pars(posterior, x)$f})
  Est_f <- rowMeans(plogis(f_sams_pf))
  plot(Sum_f$mean, Est_f, xlab = "HMC", xlim = c(0, 1), ylim = c(0, 1),
       ylab = "Pathfinder", main = paste0("est ELBO:", est_ELBO))
  abline(a = 0, b = 1)
  readline(prompt="Press [enter] to continue:")
}

