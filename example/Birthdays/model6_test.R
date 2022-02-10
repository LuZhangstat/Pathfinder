## birthday example ##
setwd("./example") # set working dir to cloned package
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(tidyverse)
library(cmdstanr)
library(posterior)
options(pillar.neg = FALSE, pillar.subtle=FALSE, pillar.sigfig=2)
library(loo)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(patchwork)
set1 <- RColorBrewer::brewer.pal(7, "Set1")
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")


#' ## Load and plot data
#' 
#' Load birthdays per day in USA 1969-1988:
data <- read_csv("./Birthdays/data/births_usa_1969.csv")

#' Add date type column for plotting
data <- data %>%
  mutate(date = as.Date("1968-12-31") + id,
         births_relative100 = births/mean(births)*100)

#' ### Plot all births
#'
#' We can see slow variation in trend, yearly pattern, and especially
#' in the later years spread to lower and higher values.
data %>%
  ggplot(aes(x=date, y=births)) + geom_point(color=set1[2]) +
  labs(x="Date", y="Relative number of births")

#' ### Plot all births as relative to mean
#'
#' To make the interpretation we switch to examine the relative
#' change, with the mean level denoted with 100.
data %>%
  ggplot(aes(x=date, y=births_relative100)) + geom_point(color=set1[2]) +
  geom_hline(yintercept=100, color='gray') +
  labs(x="Date", y="Relative births per day")

#' ### Plot mean per day of year
#'
#' We can see the generic pattern in yearly seasonal trend simply by
#' averaging over each day of year (day_of_year has numbers from 1 to
#' 366 every year with leap day being 60 and 1st March 61 also on
#' non-leap-years).
data %>%
  group_by(day_of_year2) %>%
  summarise(meanbirths=mean(births_relative100)) %>%
  ggplot(aes(x=as.Date("1986-12-31")+day_of_year2, y=meanbirths)) +
  geom_point(color=set1[2]) +
  geom_hline(yintercept=100, color='gray') +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x="Day of year", y="Relative births per day of year")

#' ### Plot mean per day of week
#' 
#' We can see the generic pattern in weekly trend simply by averaging
#' over each day of week.
data %>%
  group_by(day_of_week) %>%
  summarise(meanbirths=mean(births_relative100)) %>%
  ggplot(aes(x=day_of_week, y=meanbirths)) +
  geom_point(color=set1[2], size=4) +
  geom_hline(yintercept=100, color='gray') +
  scale_x_continuous(breaks = 1:7, labels=c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')) +
  labs(x="Day of week", y="Relative number of births of week")

## compile model and load data ##
model6 <- cmdstan_model(stan_file = "/home/luzhang/Downloads/casestudies-master/Birthdays/gpbf6.stan",
                        include_paths = "/home/luzhang/Downloads/casestudies-master/Birthdays")


standata6 <- list(x=data$id,
                  y=log(data$births_relative100),
                  N=length(data$id),
                  c_f1=1.5, # factor c of basis functions for GP for f1
                  M_f1=10, # number of basis functions for GP for f1
                  J_f2=20, # number of basis functions for periodic f2
                  day_of_week=data$day_of_week,
                  day_of_year=data$day_of_year2) # 1st March = 61 every year


## generate reference posterior samples ##
fit_stan_flag <- FALSE
if(fit_stan_flag){
  
  # use optimization to set initials #
  opt6 <- model6$optimize(data=standata6, init=0, algorithm='lbfgs',
                          history=100, tol_obj=10)
  #' Check whether parameters have reasonable values
  odraws6 <- opt6$draws()
  subset(odraws6, variable=c('intercept','sigma_','lengthscale_','sigma'), regex=TRUE)
  subset(odraws6, variable=c('beta_f3'))
  
  init6 <- sapply(c('intercept0','lengthscale_f1','lengthscale_f2',
                    'sigma_f1','sigma_f2','sigma_f4','sigma',
                    'beta_f1','beta_f2','beta_f3','beta_f4'),
                  function(variable) {as.numeric(subset(odraws6, variable=variable))})
  
  
  fit6 <- model6$sample(data=standata6,
                        chains=4, parallel_chains=4,
                        init=function() { init6 },
                        seed = 1948458383,
                        adapt_delta = 0.9,
                        thin = 100,
                        iter_warmup = 50000,
                        iter_sampling = 250000,
                        show_messages = TRUE,
                        sig_figs = 16)
  fit6$print("lp__")
  fit6$save_object(file = "./Birthdays/birthday6_ref.RDS")
  p1 <- mcmc_trace(fit6$draws("lp__"), iter1 = 1) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  #' Check whether parameters have reasonable values
  draws6 <- fit6$draws()
  summarise_draws(subset(draws6, variable=c('intercept','sigma_','lengthscale_',
                                            'sigma'), regex=TRUE))
  summarise_draws(subset(draws6, variable=c('beta_f3')))
  
  #' Compare the model to the data
  draws6 <- as_draws_matrix(draws6)
  Ef <- exp(apply(subset(draws6, variable='f'), 2, median))
  Ef1 <- apply(subset(draws6, variable='f1'), 2, median)
  Ef1 <- exp(Ef1 - mean(Ef1) + mean(log(data$births_relative100)))
  Ef2 <- apply(subset(draws6, variable='f2'), 2, median)
  Ef2 <- exp(Ef2 - mean(Ef2) + mean(log(data$births_relative100)))
  Ef_day_of_week <- apply(subset(draws6, variable='f_day_of_week'), 2, median)
  Ef_day_of_week <- exp(Ef_day_of_week - mean(Ef_day_of_week) + mean(log(data$births_relative100)))
  pf <- data %>%
    mutate(Ef = Ef) %>%
    ggplot(aes(x=date, y=births_relative100)) + geom_point(color=set1[2], alpha=0.2) +
    geom_line(aes(y=Ef), color=set1[1], alpha=0.75) +
    labs(x="Date", y="Relative number of births")
  pf1 <- data %>%
    mutate(Ef1 = Ef1) %>%
    ggplot(aes(x=date, y=births_relative100)) + geom_point(color=set1[2], alpha=0.2) +
    geom_line(aes(y=Ef1), color=set1[1]) +
    geom_hline(yintercept=100, color='gray') +
    labs(x="Date", y="Relative number of births")
  pf2 <- data %>%
    mutate(Ef2 = Ef2) %>%
    group_by(day_of_year2) %>%
    summarise(meanbirths=mean(births_relative100), meanEf2=mean(Ef2)) %>%
    ggplot(aes(x=as.Date("1987-12-31")+day_of_year2, y=meanbirths)) + geom_point(color=set1[2], alpha=0.2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    geom_line(aes(y=meanEf2), color=set1[1]) +
    geom_hline(yintercept=100, color='gray') +
    labs(x="Date", y="Relative number of births of year")
  pf3 <- ggplot(data=data, aes(x=day_of_week, y=births_relative100)) + geom_point(color=set1[2], alpha=0.2) +
    scale_x_continuous(breaks = 1:7, labels=c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')) +
    geom_line(data=data.frame(x=1:7,y=Ef_day_of_week), aes(x=x, y=Ef_day_of_week), color=set1[1]) +
    geom_hline(yintercept=100, color='gray') +
    labs(x="Date", y="Relative number of births of week")
  (pf + pf1) / (pf2 + pf3)
} 
## save the unconstrained reference posterior samples ##
model <- stan_model(file = "./Birthdays/gpbf6_rstan.stan")
posterior <- to_posterior(model, standata6)

## load reference posterior samples ##
if(fit_stan_flag){
  Stan_draws <- matrix(draws6, nrow = dim(draws6)[1])
  colnames(Stan_draws) <- dimnames(draws6)$variable
  unconstrained_draws <- unconstrain_cmd_draws(Stan_draws, posterior)
  
  lp_INV <- quantile(Stan_draws[, "lp__"], c(0.005, 0.995))
  
  save(file = "../example/Birthdays/reference_model6.RData",
       list = c("unconstrained_draws", "lp_INV"))
}else{
  load("./Birthdays/reference_model6.RData")
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
  opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, standata6,
                                  init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                  factr_tol, lmm) # plot for 8school init_bound = 15
  print(proc.time() - t)
  
  pick_samples_IR <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)
  #pick_samples_IR <- Imp_Resam_WOR(opath, n_inits = 20, seed = 1)
  pick_samples <- pick_samples_IR
  pf_fn_calls <- sapply(opath, f <- function(x){x$fn_call})
  pf_gr_calls <- sapply(opath, f <- function(x){x$gr_call})
  save(file = "./Birthdays/opath6.RData",
       list = c("opath", "pick_samples", "pf_fn_calls", "pf_gr_calls"))
}else{
  load("./Birthdays/opath6.RData")
}

## PhI ##
fit_PhI_flag = FALSE
###  random inits  ###
if(fit_PhI_flag == TRUE){
  fit_PhI <- model6$sample(
    data = standata6,
    seed = 123,
    chains = 20,
    parallel_chains = 10,
    #refresh = 0,
    save_warmup = TRUE,
    iter_warmup = 100,
    iter_sampling = 1,
    init_buffer = 100,
    term_buffer = 0,
    show_messages = TRUE,
    sig_figs = 16)
  fit_PhI$print("lp__")
  fit_PhI$save_object(file = "./Birthdays/model6_phI.RDS")
  
  p1 <- mcmc_trace(fit_PhI$draws("lp__", inc_warmup = TRUE)[10:75, ,], iter1 = 10) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  ###  record the number of iterations required to reach INV ###
  # Get the number of iterations and  leapfrogs #
  
  # Get the cost of Phase I warmup
  PhI_leapfrog_counts <- colSums(
    fit_PhI$sampler_diagnostics(inc_warmup = TRUE)[1:75, , "n_leapfrog__"])
  fit_draws_PhI <- fit_PhI$draws(inc_warmup = TRUE)
  
  last_pI_draws <- matrix(fit_draws_PhI[75, , ], 
                          nrow = dim(fit_draws_PhI)[2])
  colnames(last_pI_draws) <- dimnames(fit_draws_PhI)$variable
  unconstrained_draws_pI <- unconstrain_cmd_draws(last_pI_draws, posterior)
  PhaseI_last_draw <- unconstrained_draws_pI
  
  save(file = "./Birthdays/unconstrained_draws_model6_phI.RData",
       list = c("PhI_leapfrog_counts", "PhaseI_last_draw"))
  
} else{
  fit_PhI <- readRDS(file = "./Birthdays/model6_phI.RDS")
  load("./Birthdays/unconstrained_draws_model6_phI.RData")
}

###  inits from pf  ###
fit_PhI_pf_flag = FALSE
if(fit_PhI_pf_flag == TRUE){
  ## transform draws from pathfinder into constrained space
  init_pf_ls <- apply(pick_samples[, 1:20], 2, 
                      f <- function(x){constrain_pars(posterior, x)})
  
  fit_PhI_pf <- model6$sample(
    data = standata6,
    seed = 1234,
    init = init_pf_ls,
    chains = 20,
    parallel_chains = 10,
    #refresh = 0,
    save_warmup = TRUE,
    iter_warmup = 100,
    iter_sampling = 1,
    init_buffer = 100,
    term_buffer = 0,
    show_messages = TRUE,
    sig_figs = 16)
  fit_PhI_pf$print("lp__")
  fit_PhI_pf$save_object(file = "./Birthdays/model6_phI_pf.RDS")
  
  p1_pf <- mcmc_trace(fit_PhI_pf$draws("lp__", inc_warmup = TRUE)[10:75, ,], 
                      iter1 = 10) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1_pf)
  
  ###  record the number of iterations required to reach INV ###
  # Get the number of iterations and  leapfrogs #
  
  # Get the cost of Phase I warmup
  PhI_pf_leapfrog_counts <- colSums(
    fit_PhI_pf$sampler_diagnostics(inc_warmup = TRUE)[1:75, , "n_leapfrog__"])
  fit_draws_PhI_pf <- fit_PhI_pf$draws(inc_warmup = TRUE)
  
  last_phI_pf_draws <- matrix(fit_draws_PhI_pf[75, , ], 
                              nrow = dim(fit_draws_PhI_pf)[2])
  colnames(last_phI_pf_draws) <- dimnames(fit_draws_PhI_pf)$variable
  unconstrained_draws_pI_pf <- unconstrain_cmd_draws(last_phI_pf_draws, posterior)
  PhaseI_pf_last_draw <- unconstrained_draws_pI_pf
  
  save(file = "./Birthdays/unconstrained_draws_model6_phI_pf.RData",
       list = c("PhI_pf_leapfrog_counts", "PhaseI_pf_last_draw"))
  
} else{
  fit_PhI_pf <- readRDS(file = "./Birthdays/model6_phI_pf.RDS")
  load("./Birthdays/unconstrained_draws_model6_phI_pf.RData")
}

## check plots ##
check_dim <- c(424, 426)
for(first_check in 1:4){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  
  unconstrained_draws_PI <- PhaseI_last_draw #first_sample_draw #PhaseI_last_draw #PhaseI_last_draw[[i]] #ADVI_fullrank_draw[[i]][1:20, ]
  unconstrained_draws_PI_pf <- PhaseI_pf_last_draw
  
  plot(unconstrained_draws[, check_dim[1]], unconstrained_draws[, check_dim[2]], 
       col = alpha("grey", 0.2),
       xlim = range(unconstrained_draws[, check_dim[1]],
                    #unconstrained_draws_PI[, check_dim[1]],
                    unconstrained_draws_PI_pf[, check_dim[1]],
                    pick_samples[check_dim[1], ]
       ),
       ylim = range(unconstrained_draws[, check_dim[2]], 
                    #unconstrained_draws_PI[, check_dim[2]],
                    unconstrained_draws_PI_pf[, check_dim[2]],
                    pick_samples[check_dim[2], ]
       ),
       main = paste0("element ", check_dim[1], " and ", check_dim[2]))
  
  # points(x = unconstrained_draws_PI[, check_dim[1]],
  #        y = unconstrained_draws_PI[, check_dim[2]], col = "green",
  #        pch = 16)
  
  points(x = unconstrained_draws_PI_pf[, check_dim[1]],
         y = unconstrained_draws_PI_pf[, check_dim[2]], col = "blue",
         pch = 16)
  
  points(x = pick_samples[check_dim[1], ], y = pick_samples[check_dim[2], ],
         col = "orange", pch = 16)
  
  readline(prompt="Press [enter] to continue:")
}

## check the Wasserstein distance ##
library(transport)
ref_samples = unconstrained_draws

# pathfinder #
if(ncol(pick_samples) == 1){
  a = wpp(rbind(t(pick_samples), t(pick_samples)), 
          mass = rep(1 / 2, 2))
}else{
  a = wpp(t(pick_samples[, 1:20]), 
          mass = rep(1 / ncol(pick_samples[, 1:20]), ncol(pick_samples[, 1:20])))
}
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_pf <- wasserstein(a, b, p = 1); w_d_pf #5.45  #5.51 


# last samples of phase I #
a_phI = wpp(PhaseI_last_draw, #first_sample_draw, #PhaseI_last_draw, #[-c(2, 7, 12, 15, 16), ],
            mass = rep(1 / nrow(PhaseI_last_draw),
                       nrow(PhaseI_last_draw)))
w_d_phI <- wasserstein(a_phI, b, p = 1); w_d_phI # 9.65

# last samples of pathfinder + phase I#
a_phI_pf = wpp(PhaseI_pf_last_draw, #first_sample_draw, #PhaseI_last_draw, #[-c(2, 7, 12, 15, 16), ],
               mass = rep(1 / nrow(PhaseI_pf_last_draw),
                          nrow(PhaseI_pf_last_draw)))
w_d_phI_pf <- wasserstein(a_phI_pf, b, p = 1); w_d_phI_pf # 5.32 

## cost ##
# pf
pf_fn <- sapply(opath, f <- function(x){x$fn_call})
pf_gr <- sapply(opath, f <- function(x){x$gr_call})
sum(pf_fn);max(pf_fn)
sum(pf_gr);max(pf_gr)

#Stan phI
PhI_leapfrog_counts
PhI_pf_leapfrog_counts



## plots of fitted points ##
check_dim <- c(1, 2)  #24
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])

dta_phI <- data.frame(x = c(ref_samples[, check_dim[1]], 
                            PhaseI_last_draw[, check_dim[1]]),
                      y = c(ref_samples[, check_dim[2]],
                            PhaseI_last_draw[, check_dim[2]]),
                      label = c(rep("1", length(ref_samples[, check_dim[1]])),
                                rep("2", length(PhaseI_last_draw[, check_dim[2]]))))

dta_phI_pf <- data.frame(x = c(ref_samples[, check_dim[1]], 
                            PhaseI_pf_last_draw[, check_dim[1]]),
                      y = c(ref_samples[, check_dim[2]],
                            PhaseI_pf_last_draw[, check_dim[2]]),
                      label = c(rep("1", length(ref_samples[, check_dim[1]])),
                                rep("2", length(PhaseI_pf_last_draw[, check_dim[2]]))))

dta_pf <- data.frame(x = c(ref_samples[, check_dim[1]], 
                           pick_samples[check_dim[1], ]),
                         y = c(ref_samples[, check_dim[2]],
                               pick_samples[check_dim[2], ]),
                         label = c(rep("1", length(ref_samples[, check_dim[1]])),
                                   rep("2", length(pick_samples[check_dim[1], ]))))



p_phI <- ggplot(dta_phI, aes(x=x, y=y, color = label) ) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("1", "2"),
                     values=c(alpha("grey", 0.1), "orange"))+
  xlab("parameter 1") + ylab("parameter 2") +
  theme(legend.position='none')
p_phI

p_phI_pf <- ggplot(dta_phI_pf, aes(x=x, y=y, color = label) ) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("1", "2"),
                     values=c(alpha("grey", 0.1), "orange"))+
  xlab("parameter 1") + ylab("parameter 2") +
  theme(legend.position='none')
p_phI_pf

p_pf <- ggplot(dta_pf, aes(x=x, y=y, color = label) ) +
  geom_point(size = 3) +
  scale_color_manual(breaks = c("1", "2"),
                     values=c(alpha("grey", 0.1), "orange"))+
  xlab("parameter 1") + ylab("parameter 2") +
  theme(legend.position='none')
p_pf

ggsave("birthday_pts.eps",
       plot = p_phI,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 7.0, height = 5.0, units = "in")

ggsave("birthday_pts_pf.eps",
       plot = p_phI_pf,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 7.0, height = 5.0, units = "in")

ggsave("birthday_pts_pf2.eps",
       plot = p_pf,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 7.0, height = 5.0, units = "in")

# plot(1:length(opath[[1]]$lgnorms), opath[[1]]$lgnorms)

# check lp plots
p1_pf <- mcmc_trace(fit_PhI_pf$draws("lp__", inc_warmup = TRUE)[10:75, ,], 
                    iter1 = 10) + ylab("log-density") + xlab("iteration") +
  theme(legend.position = "none")
p1 <- mcmc_trace(fit_PhI$draws("lp__", inc_warmup = TRUE)[10:75, ,], 
                 iter1 = 10) + ylab("log-density") + xlab("iteration") +
  theme(legend.position = "none")

ggsave("birthday_p1.eps",
       plot = p1,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 7.0, height = 5.0, units = "in")

ggsave("birthday_p1_pf.eps",
       plot = p1_pf,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 7.0, height = 5.0, units = "in")


## prediction check ##
mean_brith <- mean(data$births)
f_ref_sams <- 
  apply(unconstrained_draws, 1, f <- function(x){constrain_pars(posterior, x)$f})
Est_f_ref <- apply(exp(f_ref_sams), 1, mean)
Est_f_ref <- Est_f_ref / 100 * mean_brith

f_pf_sams <- 
  apply(pick_samples[, 1:20], 2, f <- function(x){constrain_pars(posterior, x)$f})
Est_f_pf <- apply(exp(f_pf_sams), 1, mean)
Est_f_pf <- Est_f_pf / 100 * mean_brith

f_phI_sams <- 
  apply(PhaseI_last_draw, 1, f <- function(x){constrain_pars(posterior, x)$f})
Est_f_phI <- apply(exp(f_phI_sams), 1, mean)
Est_f_phI <- Est_f_phI / 100 * mean_brith


f_pf_phI_sams <- 
  apply(PhaseI_pf_last_draw, 1, f <- function(x){constrain_pars(posterior, x)$f})
#Est_f_pf_phI <- exp(apply(f_pf_phI_sams, 1, median))
Est_f_pf_phI <- apply(exp(f_pf_phI_sams), 1, mean)
Est_f_pf_phI <- Est_f_pf_phI / 100 * mean_brith

plot(Est_f_ref, Est_f_pf, xlab = "Ef from reference samples",
     ylab = "Ef from pathfinder samples")
abline(a = 0, b = 1)


plot(Est_f_ref, Est_f_phI, xlab = "Ef from reference samples",
     ylab = "Ef from last iter of Stan Phase I samples")
abline(a = 0, b = 1)

plot(Est_f_ref, Est_f_pf_phI, xlab = "Ef from reference samples",
     ylab = "Ef from last iter of Stan Phase I samples with inits from pf")
abline(a = 0, b = 1)


dat_compar <- data.frame(Est_f_ref = Est_f_ref, Est_f_pf = Est_f_pf)
p_compar <- ggplot(dat_compar, aes(x=Est_f_ref, y=Est_f_pf)) + 
  geom_point(alpha = 0.2) +
  xlim(6000, 13000) + ylim(6000, 13000) + 
  geom_line(data = data.frame(x = c(6000, 13000), y = c(6000, 13000)),
            aes(x = x, y = y), color = "red", linetype="dashed") +
  labs(x="number of births (reference)", y="number of births (Pathfinder)") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
p_compar
ggsave("birth_pred_compar_pf.eps", #"8-school_opt_tr22.eps"
       plot = p_compar,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 4.0, height = 4.0, units = "in")


dat_compar <- data.frame(Est_f_ref = Est_f_ref, Est_f_phI = Est_f_phI)
p_compar2 <- ggplot(dat_compar, aes(x=Est_f_ref, y=Est_f_phI)) + 
  geom_point(alpha = 0.2) +
  xlim(6000, 13000) + ylim(6000, 13000) + 
  geom_line(data = data.frame(x = c(6000, 13000), y = c(6000, 13000)),
            aes(x = x, y = y), color = "red", linetype="dashed") +
  labs(x="number of births (reference)", y="number of births (NUTS)") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
p_compar2
ggsave("birth_pred_compar_phI.eps", #"8-school_opt_tr22.eps"
       plot = p_compar2,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 4.0, height = 4.0, units = "in")

dat_compar <- data.frame(Est_f_ref = Est_f_ref, Est_f_pf_phI = Est_f_pf_phI)
p_compar3 <- ggplot(dat_compar, aes(x=Est_f_ref, y=Est_f_pf_phI)) + 
  geom_point(alpha = 0.2) +
  xlim(6000, 13000) + ylim(6000, 13000) + 
  geom_line(data = data.frame(x = c(6000, 13000), y = c(6000, 13000)),
            aes(x = x, y = y), color = "red", linetype="dashed") +
  labs(x="number of births (reference)", y="number of births (pathfinder + NUTS)") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
p_compar3
ggsave("birth_pred_compar_pf_phI.eps", #"8-school_opt_tr22.eps"
       plot = p_compar3,
       device = cairo_ps,
       path = "./Birthdays/",
       width = 4.0, height = 4.0, units = "in")


mean(PhI_leapfrog_counts)
# 13488.05
mean(pf_gr_calls)
# 1286.9
mean(pf_fn_calls)
# 6386.75
mean(PhI_leapfrog_counts) / mean(pf_gr_calls + (pf_fn_calls - pf_gr_calls  + 100)/4)
# 5.2


