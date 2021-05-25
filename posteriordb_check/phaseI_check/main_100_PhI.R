rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
L = 75
M = 100
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234

load(file = "../results/lp_posteriordb_LBFGS_h6.RData")

# preallocate results #
PhaseI_last_draw <- list() # Get the last samples of Phase I
PhI_leapfrog_counts <- array(data = NA, dim = c(M, length(model_record)))

i = which(model_record == 27)
for(i in 49:49){ #20 length(model_record)
  gc()
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  ### setup initials ###
  data = get_data(po)
  posterior <- to_posterior(model, data)

  ###  run Stan with a long Phase I warmup time  ###
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ###  random inits  ###
  suppressWarnings(
    fit <- mod$sample(
      data = data,
      seed = 12345, #123 for i = 21
      chains = M,
      parallel_chains = 5,
      refresh = 0,
      save_warmup = TRUE,
      iter_warmup = L,
      iter_sampling = 0,
      init_buffer = L,
      term_buffer = 0,
      show_messages = FALSE,
      sig_figs = 18
    ))
  p1 <- mcmc_trace(fit$draws("lp__", inc_warmup = TRUE)[60:75, ,], iter1 = 60) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  ## record the initial and the last sample of phase I ##
  fit_draws <- fit$draws(inc_warmup = TRUE)
  last_pI_draws <- matrix(fit_draws[75, , ], 
                          nrow = dim(fit_draws)[2])
  colnames(last_pI_draws) <- dimnames(fit_draws)$variable
  unconstrained_draws <- unconstrain_cmd_draws(last_pI_draws, posterior)
  PhaseI_last_draw[[i]] <- unconstrained_draws # record the last sample of phase I
  
  # Get the cost of Phase I warmup
  PhI_leapfrog_counts[, i] <- colSums(
    fit$sampler_diagnostics(inc_warmup = TRUE)[1:75, , "n_leapfrog__"])
  
}

save(file = "../results/PhI_100_h10.RData",
     list = c("PhaseI_last_draw", "PhI_leapfrog_counts"))

load("../results/PhI_100_h10.RData")
length(PhaseI_last_draw)


## old code ##
# 10 random samples #
MC = 100
load("../results/PhI_100_h10.RData")
set.seed(123)

phI_10_draws <- list()

which(model_record == 48)
t_0 <- proc.time()
for(i in 1:49){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  # pick model
  
  phI_10_draws[[i]] <- list()
  for(j in 1:MC){
    set.seed(j)
    pick_pf <- sample.int(MC, 10, replace = FALSE)
    phI_10_draws[[i]][[j]] <- PhaseI_last_draw[[i]][pick_pf, ]
  }
}
proc.time() - t_0

save(file = "../results/phI_10_draws.RData",
     list = c("phI_10_draws"))


